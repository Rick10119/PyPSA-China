# scripts/build_co2_totals.py
import pandas as pd

def build_co2_totals():
    """
    Create CO2 totals reference data for China's electricity and heating sectors.
    Values are based on 2023 emissions data.
    """
    
    # Create DataFrame with CO2 emissions data
    # Units: tCO2/year
    co2_data = pd.DataFrame({
        'co2': [
            5.917263355e9  # Total emissions from electricity (5.288987673e9) and heating (0.628275682e9) sectors, 与论文中一致
        ]
    })
    
    # Save to HDF5 file
    co2_data.to_hdf(snakemake.output[0], key='co2')

if __name__ == "__main__":
    build_co2_totals()