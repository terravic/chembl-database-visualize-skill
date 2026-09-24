# SPDX-License-Identifier: Apache-2.0
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Test runner executing unit, visualization, and pipeline integration tests."""

from __future__ import annotations

import os
import sys
import unittest

if __name__ == "__main__":
  test_dir = os.path.dirname(os.path.abspath(__file__))
  suite = unittest.defaultTestLoader.discover(start_dir=test_dir, pattern="test_*.py")
  runner = unittest.TextTestRunner(verbosity=2)
  result = runner.run(suite)
  sys.exit(0 if result.wasSuccessful() else 1)
