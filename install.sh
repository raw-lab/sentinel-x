#!/usr/bin/env bash
# Build Sentinel-X (release) and install the `sentinel-x` binary onto your PATH.
#
#   ./install.sh            # install to ~/.local/bin (no sudo)
#   PREFIX=/usr/local ./install.sh   # system-wide (may need sudo)
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$here"

echo "[sentinel-x] building release binary..."
# single crate (sentinel-x); --locked uses the tested dependency versions,
# which build with rustc >= 1.75 (e.g. Ubuntu's apt cargo)
"${CARGO:-cargo}" build --release --locked

bin="$here/target/release/sentinel-x"
[ -x "$bin" ] || { echo "build did not produce $bin" >&2; exit 1; }

if [ -n "${PREFIX:-}" ]; then
    dest="$PREFIX/bin"
else
    dest="$HOME/.local/bin"
fi
mkdir -p "$dest"
install -m 0755 "$bin" "$dest/sentinel-x"
echo "[sentinel-x] installed to $dest/sentinel-x"

case ":$PATH:" in
    *":$dest:"*) echo "[sentinel-x] $dest is on your PATH — run: sentinel-x help" ;;
    *) echo "[sentinel-x] NOTE: add $dest to your PATH, e.g.:"
       echo "           echo 'export PATH=\"$dest:\$PATH\"' >> ~/.bashrc && source ~/.bashrc" ;;
esac
