import matplotlib.pyplot as plt
import numpy as np

# Data
time_durations = ['24 hours', '1 week', '2 weeks', '1 month']
storage_percentages = [20, 67, 90, 100]

# Create figure and axis
plt.figure(figsize=(10, 6))
bars = plt.bar(time_durations, storage_percentages, color='#1f77b4')

# Add horizontal line at 67%
plt.axhline(y=67, color='red', linestyle='--', linewidth=2)

# Customize the plot
plt.title('Storage Space Requirements for Aluminum Products', fontsize=14, pad=20)
plt.xlabel('Maximum Storage Duration', fontsize=12)
plt.ylabel('Reduced Energy Cost (%)', fontsize=12)

# Add value labels on top of each bar
for bar in bars:
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2., height,
             f'{int(height)}%',
             ha='center', va='bottom')

plt.show()
# Adjust layout and save
plt.tight_layout()
# plt.savefig('examples/plots/storage_requirements.png', dpi=300, bbox_inches='tight')
# plt.close() 
