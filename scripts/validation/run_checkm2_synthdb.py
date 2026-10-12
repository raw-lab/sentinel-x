"""Run CheckM2 1.1.0 unmodified except for the DIAMOND-database checksum
check, which would reject a synthetic test database."""
import sys
sys.path.insert(0, '/home/claude/refsrc/checkm2')
from checkm2 import versionControl
versionControl.VersionControl.checksum_version_validate_DIAMOND = lambda self, location=None: True
from checkm2.main import main
sys.argv = ['checkm2'] + sys.argv[1:]
main()
