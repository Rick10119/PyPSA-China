import matplotlib.pyplot as plt
import numpy as np

# Data
p_min_pu = ['0.1', '0.5', '0.9']
reduced_energy_cost = np.array([6.7, 6.0, 5.05])/91.1 * 100
start_up_cost = np.array([0.0, 1.2, 1.2])/91.1 * 100

# Set up the figure and axis
plt.figure(figsize=(10, 6))
x = np.arange(len(p_min_pu))
width = 0.35

# Create bars
bars1 = plt.bar(x - width/2, reduced_energy_cost, width, label='Reduced Energy Cost', color='#1f77b4')
bars2 = plt.bar(x + width/2, start_up_cost, width, label='Start-up Cost', color='#ff7f0e')

# Customize the plot
plt.title('Flexibility Impact on Energy and Start-up Costs', fontsize=14, pad=20)
plt.xlabel('Minimum Power Consumption (p.u.)', fontsize=12)
plt.ylabel('Cost Reduction (%)', fontsize=12)
plt.xticks(x, p_min_pu)
plt.legend()

# Add value labels on top of each bar
for bars in [bars1, bars2]:
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height,
                f'{height:.1f}%',
                ha='center', va='bottom')

# Adjust layout and save
plt.show()
plt.tight_layout()
plt.savefig('examples/plots/flexibility_impact.png', dpi=300, bbox_inches='tight')
plt.close() 
