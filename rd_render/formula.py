"""
Provides validation for chemical formulas:
- Uses regex to check element symbols and counts (e.g., C9H8O4, C2H5OH)
- Ensures element counts are reasonable (≤ 100 atoms per element)
- Returns True if the formula is syntactically valid and within limits, otherwise False
"""

import re

def is_valid_formula(formula):
    # Basic chemical formula validator: matches strings like C9H8O4 or C2H5OH
    pattern = r'^([A-Z][a-z]?\d*)+$'
    if not re.match(pattern, formula):
        return False
    # Check that element counts are reasonable (no crazy numbers)
    tokens = re.findall(r'[A-Z][a-z]?(\d*)', formula)
    for num in tokens:
        if num:
            try:
                if int(num) > 100:
                    return False
            except ValueError:
                return False
    return True
