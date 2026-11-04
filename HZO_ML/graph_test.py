from periodic_colors import get_color

c = get_color("Fe")  # Example usage
print(f"The color for Iron (Fe) is: {c}")

for i in range(1, 119):
    symbol = get_color(i)
    print(f"Element {i} color: {symbol}")
