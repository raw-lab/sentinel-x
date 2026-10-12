//! # sentinel-gbm
//!
//! A pure-Rust loader and evaluator for **LightGBM text models** — the exact
//! `.gbm` files CheckM2 ships (`general_model_COMP.gbm`, `model_CONT.gbm`).
//! It parses the model into flat decision-tree arrays and reproduces LightGBM's
//! prediction without any Python, C LightGBM, or FFI, so CheckM2-style
//! completeness/contamination inference runs entirely inside Sentinel.
//!
//! The supported model dialect is the one CheckM2 uses: dense numerical
//! features, `decision_type = 2` (numerical `<=` split, default-left on
//! missing), single output per iteration, and the `regression` /
//! `regression_sqrt` objectives.
//!
//! ```no_run
//! use crate::gbm::Booster;
//! # fn demo(model_text: &str, features: &[f64]) -> anyhow::Result<()> {
//! let booster = Booster::parse(model_text)?;
//! let yhat = booster.predict(features);
//! # let _ = yhat; Ok(()) }
//! ```

// Evaluator API: batch helpers and model metadata are used by tests only.
#![allow(dead_code)]

use anyhow::{anyhow, bail, Context, Result};
use std::path::Path;

/// The boosting objective, controlling the final output transform.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Objective {
    /// Plain L2 regression: output is the raw sum of leaf values.
    Regression,
    /// `regression_sqrt`: the model is trained in signed-sqrt space; the final
    /// output is `sign(raw) * raw^2`.
    RegressionSqrt,
}

/// One decision tree, stored as parallel flat arrays (LightGBM layout).
///
/// Internal node `i` tests `feature[ split_feature[i] ] <= threshold[i]`,
/// branching to `left_child[i]` (true) or `right_child[i]` (false). A child
/// `c < 0` denotes leaf index `-c - 1`.
#[derive(Clone, Debug)]
pub struct Tree {
    pub split_feature: Vec<usize>,
    pub threshold: Vec<f64>,
    pub left_child: Vec<i32>,
    pub right_child: Vec<i32>,
    pub leaf_value: Vec<f64>,
    /// `decision_type[i]` bitmask; bit 1 (`& 2`) = default-left on missing.
    pub decision_type: Vec<u8>,
}

impl Tree {
    /// Evaluate this tree for a dense feature vector, returning the leaf value.
    #[inline]
    pub fn predict(&self, features: &[f64]) -> f64 {
        if self.left_child.is_empty() {
            return self.leaf_value.first().copied().unwrap_or(0.0);
        }
        let mut node: i32 = 0;
        loop {
            let i = node as usize;
            let f = self.split_feature[i];
            let mut fval = features.get(f).copied().unwrap_or(0.0);
            // LightGBM's NumericalDecision: missing type in bits 2-3
            // (0 none, 1 zero, 2 NaN), default direction in bit 1
            let dt = self.decision_type[i];
            let missing = (dt >> 2) & 3;
            if fval.is_nan() && missing != 2 {
                fval = 0.0;
            }
            let go_left = if (missing == 1 && fval.abs() <= 1e-35) || (missing == 2 && fval.is_nan()) {
                dt & 2 != 0
            } else {
                fval <= self.threshold[i]
            };
            node = if go_left { self.left_child[i] } else { self.right_child[i] };
            if node < 0 {
                let leaf = (-node - 1) as usize;
                return self.leaf_value[leaf];
            }
        }
    }
}

/// A full gradient-boosted ensemble loaded from a LightGBM text model.
#[derive(Clone, Debug)]
pub struct Booster {
    pub trees: Vec<Tree>,
    pub objective: Objective,
    pub max_feature_idx: usize,
    pub num_class: usize,
}

impl Booster {
    /// Load and parse a model from a file path.
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        let text = std::fs::read_to_string(path.as_ref())
            .with_context(|| format!("reading LightGBM model {}", path.as_ref().display()))?;
        Self::parse(&text)
    }

    /// Parse a model from an in-memory LightGBM text dump.
    pub fn parse(text: &str) -> Result<Self> {
        let mut objective = Objective::Regression;
        let mut max_feature_idx = 0usize;
        let mut num_class = 1usize;

        let mut blocks: Vec<Vec<&str>> = Vec::new();
        let mut cur: Vec<&str> = Vec::new();
        let mut in_tree = false;

        for line in text.lines() {
            if line.starts_with("Tree=") {
                if !cur.is_empty() {
                    blocks.push(std::mem::take(&mut cur));
                }
                in_tree = true;
                cur.push(line);
                continue;
            }
            if in_tree {
                if line.trim().is_empty() {
                    continue;
                }
                if line.starts_with("end of trees") || line.starts_with("feature_importances:") || line.starts_with("parameters:")
                {
                    in_tree = false;
                    continue;
                }
                cur.push(line);
            } else if let Some(v) = line.strip_prefix("objective=") {
                let head = v.split_whitespace().next().unwrap_or("");
                objective = if head == "regression_sqrt" || (head == "regression" && v.contains("sqrt")) {
                    Objective::RegressionSqrt
                } else if head == "regression" || head == "regression_l2" || head == "l2" {
                    Objective::Regression
                } else {
                    bail!("unsupported LightGBM objective '{v}'");
                };
            } else if let Some(v) = line.strip_prefix("max_feature_idx=") {
                max_feature_idx = v.trim().parse().unwrap_or(0);
            } else if let Some(v) = line.strip_prefix("num_class=") {
                num_class = v.trim().parse().unwrap_or(1);
            }
        }
        if !cur.is_empty() {
            blocks.push(cur);
        }
        if blocks.is_empty() {
            bail!("no Tree= blocks found in LightGBM model");
        }
        if num_class != 1 {
            bail!("multi-class LightGBM models are not supported");
        }

        let mut trees = Vec::with_capacity(blocks.len());
        for blk in &blocks {
            trees.push(parse_tree(blk)?);
        }

        Ok(Booster { trees, objective, max_feature_idx, num_class })
    }

    /// Raw ensemble sum (sum of leaf values across all trees).
    #[inline]
    pub fn raw_score(&self, features: &[f64]) -> f64 {
        self.trees.iter().map(|t| t.predict(features)).sum()
    }

    /// Number of input features the model expects.
    pub fn n_features(&self) -> usize {
        self.max_feature_idx + 1
    }

    /// Final prediction, applying the objective transform.
    #[inline]
    pub fn predict(&self, features: &[f64]) -> f64 {
        debug_assert!(features.len() >= self.n_features(), "feature vector shorter than the model's");
        let raw = self.raw_score(features);
        match self.objective {
            Objective::Regression => raw,
            Objective::RegressionSqrt => raw.signum() * raw * raw,
        }
    }

    /// Predict for a batch of feature vectors.
    pub fn predict_batch(&self, rows: &[Vec<f64>]) -> Vec<f64> {
        rows.iter().map(|r| self.predict(r)).collect()
    }

    /// Number of boosting iterations (trees).
    pub fn n_trees(&self) -> usize {
        self.trees.len()
    }
}

fn parse_usize_list(line: &str, key: &str) -> Result<Vec<usize>> {
    let v = line.strip_prefix(key).ok_or_else(|| anyhow!("expected '{key}' line"))?;
    v.split_whitespace().map(|t| t.parse::<usize>().map_err(|e| anyhow!("{key}: {e}"))).collect()
}
fn parse_i32_list(line: &str, key: &str) -> Result<Vec<i32>> {
    let v = line.strip_prefix(key).ok_or_else(|| anyhow!("expected '{key}' line"))?;
    v.split_whitespace().map(|t| t.parse::<i32>().map_err(|e| anyhow!("{key}: {e}"))).collect()
}
fn parse_f64_list(line: &str, key: &str) -> Result<Vec<f64>> {
    let v = line.strip_prefix(key).ok_or_else(|| anyhow!("expected '{key}' line"))?;
    v.split_whitespace().map(|t| t.parse::<f64>().map_err(|e| anyhow!("{key}: {e}"))).collect()
}
fn parse_u8_list(line: &str, key: &str) -> Result<Vec<u8>> {
    let v = line.strip_prefix(key).ok_or_else(|| anyhow!("expected '{key}' line"))?;
    v.split_whitespace().map(|t| t.parse::<u8>().map_err(|e| anyhow!("{key}: {e}"))).collect()
}

fn parse_tree(block: &[&str]) -> Result<Tree> {
    let find = |key: &str| block.iter().find(|l| l.starts_with(key)).copied();

    let num_leaves: usize = find("num_leaves=")
        .and_then(|l| l.strip_prefix("num_leaves=").and_then(|v| v.trim().parse().ok()))
        .ok_or_else(|| anyhow!("tree missing num_leaves"))?;

    if num_leaves <= 1 {
        let leaf_value = find("leaf_value=").map(|l| parse_f64_list(l, "leaf_value=")).transpose()?.unwrap_or_else(|| vec![0.0]);
        return Ok(Tree {
            split_feature: vec![],
            threshold: vec![],
            left_child: vec![],
            right_child: vec![],
            leaf_value,
            decision_type: vec![],
        });
    }

    let split_feature =
        parse_usize_list(find("split_feature=").ok_or_else(|| anyhow!("missing split_feature"))?, "split_feature=")?;
    let threshold = parse_f64_list(find("threshold=").ok_or_else(|| anyhow!("missing threshold"))?, "threshold=")?;
    let left_child = parse_i32_list(find("left_child=").ok_or_else(|| anyhow!("missing left_child"))?, "left_child=")?;
    let right_child = parse_i32_list(find("right_child=").ok_or_else(|| anyhow!("missing right_child"))?, "right_child=")?;
    let leaf_value = parse_f64_list(find("leaf_value=").ok_or_else(|| anyhow!("missing leaf_value"))?, "leaf_value=")?;
    let decision_type = match find("decision_type=") {
        Some(l) => parse_u8_list(l, "decision_type=")?,
        None => vec![2u8; split_feature.len()],
    };

    let n = split_feature.len();
    if threshold.len() != n || left_child.len() != n || right_child.len() != n || decision_type.len() != n {
        bail!("tree internal-node arrays have inconsistent lengths");
    }
    if leaf_value.len() != num_leaves || n + 1 != num_leaves {
        bail!("tree has {} leaves but {} leaf values / {} splits", num_leaves, leaf_value.len(), n);
    }
    if decision_type.iter().any(|d| d & 1 != 0) {
        bail!("categorical splits are not supported");
    }
    for &c in left_child.iter().chain(right_child.iter()) {
        let ok = if c < 0 { ((-c - 1) as usize) < num_leaves } else { (c as usize) < n && c > 0 };
        if !ok {
            bail!("tree child index {c} out of range");
        }
    }

    Ok(Tree { split_feature, threshold, left_child, right_child, leaf_value, decision_type })
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOY: &str = "\
tree
version=v3
num_class=1
num_tree_per_iteration=1
max_feature_idx=1
objective=regression

Tree=0
num_leaves=2
num_cat=0
split_feature=0
threshold=0.5
decision_type=2
left_child=-1
right_child=-2
leaf_value=10 20
shrinkage=1


Tree=1
num_leaves=2
num_cat=0
split_feature=1
threshold=1.5
decision_type=2
left_child=-1
right_child=-2
leaf_value=1 2
shrinkage=0.1
";

    #[test]
    fn parses_and_predicts() {
        let b = Booster::parse(TOY).unwrap();
        assert_eq!(b.n_trees(), 2);
        assert_eq!(b.objective, Objective::Regression);
        assert!((b.predict(&[0.0, 0.0]) - 11.0).abs() < 1e-9);
        assert!((b.predict(&[1.0, 2.0]) - 22.0).abs() < 1e-9);
        assert!((b.predict(&[1.0, 0.0]) - 21.0).abs() < 1e-9);
    }

    #[test]
    fn sqrt_objective_transforms() {
        let t = TOY.replace("objective=regression", "objective=regression_sqrt");
        let b = Booster::parse(&t).unwrap();
        assert!((b.predict(&[0.0, 0.0]) - 121.0).abs() < 1e-6);
    }
}
