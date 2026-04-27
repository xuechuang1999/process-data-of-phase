# periodic_color.py
# Jmol CPK 颜色标准 - 包含元素 1 (H) - 118 (Og)
# 适用于 Matplotlib, PyMatGen, ASE 等

element_colors = {
    # Period 1
    "H":  "#FFFFFF", "He": "#D9FFFF",
    
    # Period 2
    "Li": "#CC80FF", "Be": "#C2FF00", "B":  "#FFB5B5", "C":  "#909090", 
    "N":  "#3050F8", "O":  "#FF0D0D", "F":  "#90E050", "Ne": "#B3E3F5",
    
    # Period 3
    "Na": "#AB5CF2", "Mg": "#8AFF00", "Al": "#BFA6A6", "Si": "#F0C8A0", 
    "P":  "#FF8000", "S":  "#FFFF30", "Cl": "#1FF01F", "Ar": "#80D1E3",
    
    # Period 4
    "K":  "#8F40D4", "Ca": "#3DFF00", "Sc": "#E6E6E6", "Ti": "#BFC2C7", 
    "V":  "#A6A6AB", "Cr": "#8A99C7", "Mn": "#9C7AC7", "Fe": "#E06633", 
    "Co": "#F090A0", "Ni": "#50D050", "Cu": "#C88033", "Zn": "#7D80B0", 
    "Ga": "#C28F8F", "Ge": "#668F8F", "As": "#BD80E3", "Se": "#FFA100", 
    "Br": "#A62929", "Kr": "#5CB8D1",
    
    # Period 5
    "Rb": "#702EB0", "Sr": "#00FF00", "Y":  "#94FFFF", "Zr": "#94E0E0", 
    "Nb": "#73C2C9", "Mo": "#54B5B5", "Tc": "#3B9E9E", "Ru": "#248F8F", 
    "Rh": "#0A7D8C", "Pd": "#006985", "Ag": "#C0C0C0", "Cd": "#FFD98F", 
    "In": "#A67573", "Sn": "#668080", "Sb": "#9E63B5", "Te": "#D47A00", 
    "I":  "#940094", "Xe": "#429EB0",
    
    # Period 6 (including Lanthanides)
    "Cs": "#57178F", "Ba": "#00C900", "La": "#70D4FF", "Ce": "#FFFFC7", 
    "Pr": "#D9FFC7", "Nd": "#C7FFC7", "Pm": "#A3FFC7", "Sm": "#8FFFC7", 
    "Eu": "#61FFC7", "Gd": "#45FFC7", "Tb": "#30FFC7", "Dy": "#1FFFC7", 
    "Ho": "#00FF9C", "Er": "#00E675", "Tm": "#00D452", "Yb": "#00BF38", 
    "Lu": "#00AB24", "Hf": "#4DC2FF", "Ta": "#4DA6FF", "W":  "#2194D6", 
    "Re": "#267DAB", "Os": "#266696", "Ir": "#175487", "Pt": "#D0D0E0", 
    "Au": "#FFD123", "Hg": "#B8B8D0", "Tl": "#A6544D", "Pb": "#575961", 
    "Bi": "#9E4FB5", "Po": "#AB5C00", "At": "#754F45", "Rn": "#428296",
    
    # Period 7 (including Actinides)
    "Fr": "#420066", "Ra": "#007D00", "Ac": "#70ABFA", "Th": "#00BAFF", 
    "Pa": "#00A1FF", "U":  "#008FFF", "Np": "#0080FF", "Pu": "#006BFF", 
    "Am": "#545CF2", "Cm": "#785CE3", "Bk": "#8A4FE3", "Cf": "#A136D4", 
    "Es": "#B31FD4", "Fm": "#B31FBA", "Md": "#B30DA6", "No": "#BD0D87", 
    "Lr": "#C70066", "Rf": "#CC0059", "Db": "#D1004F", "Sg": "#D90045", 
    "Bh": "#E00038", "Hs": "#E6002E", "Mt": "#EB0026", 
    
    # Trans-actinides (Super heavy elements) - Jmol Extended
    "Ds": "#ED0023", "Rg": "#F00021", "Cn": "#F2001F", "Nh": "#F4001C", 
    "Fl": "#F6001A", "Mc": "#F80017", "Lv": "#FA0014", "Ts": "#FC0012", 
    "Og": "#FC000F" 
}

# 辅助函数：通过原子序数获取颜色（可选）
# 如果你的 VASP POTCAR 只给出了 Z 值，这个会很有用
z_to_symbol = {
    1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
    11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar',
    19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
    31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr',
    37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd',
    49: 'In', 50: 'Sn', 51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe',
    55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu',
    72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg',
    81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn',
    87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th', 91: 'Pa', 92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es', 100: 'Fm', 101: 'Md', 102: 'No', 103: 'Lr',
    104: 'Rf', 105: 'Db', 106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt',
    110: 'Ds', 111: 'Rg', 112: 'Cn', 113: 'Nh', 114: 'Fl', 115: 'Mc', 116: 'Lv', 117: 'Ts', 118: 'Og'
}

def get_color(identifier):
    """
    输入: 元素符号(str) 如 'Fe' 或 原子序数(int) 如 26
    输出: Hex颜色代码
    """
    if isinstance(identifier, int):
        symbol = z_to_symbol.get(identifier)
        return element_colors.get(symbol, "#808080")
    else:
        return element_colors.get(identifier, "#808080")