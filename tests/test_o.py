#!/usr/bin/env python3
"""test_o.py -- OLCAO test suite entry point.

Run from the repository root:

    python src/tests/test_o.py               # all tests
    python src/tests/test_o.py -m unit       # unit tests only
    python src/tests/test_o.py -m integration
    python src/tests/test_o.py -v            # verbose
    python src/tests/test_o.py -v -m unit    # verbose unit tests

All arguments are forwarded directly to pytest.

If OLCAO_DATA is not set, the script looks for <repo_root>/share/ and
sets the variable automatically when that directory exists.
"""

import os
import sys


def main():
    import pytest

    tests_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(os.path.dirname(tests_dir))

    # Auto-detect OLCAO_DATA when the caller has not set it.
    if not os.environ.get('OLCAO_DATA'):
        candidate = os.path.join(repo_root, 'share')
        if os.path.isdir(candidate):
            os.environ['OLCAO_DATA'] = candidate
            print(f'[test_o] OLCAO_DATA not set; using {candidate}')
        else:
            print(
                '[test_o] Warning: OLCAO_DATA is not set and '
                f'{candidate} was not found.  Tests that need the '
                'element database will be skipped.',
                file=sys.stderr,
            )

    # If the caller supplies no arguments, run the whole tests directory
    # (pytest will use testpaths from pyproject.toml).
    args = sys.argv[1:] if len(sys.argv) > 1 else []
    sys.exit(pytest.main(args))


if __name__ == '__main__':
    main()
