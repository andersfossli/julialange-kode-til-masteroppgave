
import pandas as pd
import numpy as np

def standardize_reactor_type(reactor_type):
    """
    Standardize reactor type names to basic technology categories
    
    Parameters:
    -----------
    reactor_type : str
        Original reactor type from data
        
    Returns:
    --------
    str : Standardized reactor type
    """
    # Convert to uppercase and strip whitespace
    rtype = str(reactor_type).upper().strip()
    
    # PWR variants
    pwr_types = [
        'PWR', 'EPR', 'EPR (PWR)', 'EPR-1750', 
        'AP1000', 'APR1400', 
        'VVER1200 (PWR)', 'VVER', 'VVER1200',
        'IPWR'  # Integral PWR
    ]
    
    # BWR variants
    bwr_types = ['BWR']
    
    # HTR variants
    htr_types = ['HTR', 'HTR/GFR', 'HTGR']
    
    # SFR variants
    sfr_types = ['SFR', 'BN-800 (SFR)', 'BN-800', 'BN800']
    
    # MSR variants
    msr_types = ['MSR', 'MSFR']
    
    # LFR variants
    lfr_types = ['LFR']
    
    # MR variants
    mr_types = ['MR']
    
    # Check which category it belongs to
    if any(pwr in rtype for pwr in pwr_types):
        return 'PWR'
    elif any(bwr in rtype for bwr in bwr_types):
        return 'BWR'
    elif any(htr in rtype for htr in htr_types):
        return 'HTR'
    elif any(sfr in rtype for sfr in sfr_types):
        return 'SFR'
    elif any(msr in rtype for msr in msr_types):
        return 'MSR'
    elif any(lfr in rtype for lfr in lfr_types):
        return 'LFR'
    elif any(mr in rtype for mr in mr_types):
        return 'MR'
    else:
        print(f"⚠ Warning: Unknown reactor type '{reactor_type}', keeping as-is")
        return reactor_type


def assign_tech_category(row):
    """
    Assign combined technology category (Size-Type)
    
    Parameters:
    -----------
    row : pd.Series
        Row from dataframe with 'size_category' and 'reactor_type'
    
    Returns:
    --------
    str : Combined category like 'SMR-PWR', 'Large-BWR', etc.
    """
    size = row.get('size_category', 'Unknown')
    rtype = row.get('reactor_type', 'Unknown')
    
    return f"{size}-{rtype}"


def extract_reactor_data(excel_file, output_csv='reactor_data.csv'):
    """
    Extract reactor data from Excel for LCOE Monte Carlo simulation
    
    Parameters:
    -----------
    excel_file : str
        Path to your Excel file
    output_csv : str
        Output CSV filename
    """
    
    # Read the specific sheet
    df = pd.read_excel(excel_file, sheet_name='Datasheet', skiprows=2)
    
    print("Available columns in Excel:")
    for i, col in enumerate(df.columns, 1):
        print(f"{i}. '{col}'")
    print("\n")
    
    # Define column mapping for LCOE simulation
    column_mapping = {
        # Identifiers
        'Country': 'country',
        'Project': 'name',
        'Type': 'reactor_type',
        'FOAK/NOAK*': 'foak_noak',
        'Large medium micro': 'size_category',
        
        # Core economic parameters
        'Capacity (net MWe)': 'capacity_mwe',
        'OCC (USD2020/kW)': 'occ_usd_per_kw',
        
        # Construction time parameters
        'planned Construction time (y)': 'construction_time_planned',
        'Actual / total build time (y)': 'construction_time_actual',
        
        # Operational parameters
        'Lifetime (y)': 'lifetime_years',
        'Capacity factor (%)': 'capacity_factor_pct',
        
        # Cost parameters
        'OPEX fixed (USD2020/MW-yr)': 'opex_fixed_usd_per_mw_yr',
        'OPEX variable (USD2020/MWh)': 'opex_variable_usd_per_mwh',
        'Fuel (USD2020/MWh)': 'fuel_usd_per_mwh',
        'Waste (USD2020/MWh)': 'waste_usd_per_mwh',
        'Decom (% of CAPEX)': 'decom_pct_of_capex',
        
        # Optional but useful
        'WACC (real %)': 'wacc_pct',
        'Year': 'data_year'
    }
    
    # Check which columns exist
    available_cols = {k: v for k, v in column_mapping.items() if k in df.columns}
    missing_cols = set(column_mapping.keys()) - set(available_cols.keys())
    
    if missing_cols:
        print(f"Warning: These columns not found: {missing_cols}\n")
    
    if not available_cols:
        print("ERROR: No matching columns found!")
        return None
    
    # Select and rename columns
    df_output = df[list(available_cols.keys())].copy()
    df_output = df_output.rename(columns=available_cols)
    
    # Clean numeric columns
    numeric_cols = [
        'capacity_mwe', 'occ_usd_per_kw', 
        'construction_time_planned', 'construction_time_actual',
        'lifetime_years', 'capacity_factor_pct',
        'opex_fixed_usd_per_mw_yr', 'opex_variable_usd_per_mwh',
        'fuel_usd_per_mwh', 'waste_usd_per_mwh', 
        'decom_pct_of_capex', 'wacc_pct'
    ]
    
    for col in numeric_cols:
        if col in df_output.columns:
            # Remove spaces, dashes, and convert to numeric
            df_output[col] = pd.to_numeric(
                df_output[col].astype(str)
                    .str.replace(' ', '')
                    .str.replace('–', '')
                    .str.replace('-', '')
                    .str.replace('x', ''),  # Remove 'x' markers
                errors='coerce'
            )
    
    # === NEW: Standardize reactor types ===
    print("--- Reactor Type Standardization ---")
    if 'reactor_type' in df_output.columns:
        print("Original reactor types:")
        print(df_output['reactor_type'].value_counts())
        
        df_output['reactor_type'] = df_output['reactor_type'].apply(standardize_reactor_type)
        
        print("\nStandardized reactor types:")
        print(df_output['reactor_type'].value_counts())
        print()
    
    # === NEW: Assign combined tech category ===
    if 'size_category' in df_output.columns and 'reactor_type' in df_output.columns:
        df_output['tech_category'] = df_output.apply(assign_tech_category, axis=1)
        print("--- Technology Categories Created ---")
        print(df_output['tech_category'].value_counts())
        print()
    
    # Data validation and cleaning
    print("--- Data Cleaning ---")
    
    # Handle capacity factor ranges (e.g., "85-90" -> use midpoint)
    if 'capacity_factor_pct' in df_output.columns:
        # For now, if it's a range string, we'll set it to NaN and handle separately
        # You might want to parse ranges like "85-90" -> 87.5
        pass
    
    # Calculate absolute OCC (total capital cost)
    if 'occ_usd_per_kw' in df_output.columns and 'capacity_mwe' in df_output.columns:
        df_output['occ_total_million_usd'] = (
            df_output['occ_usd_per_kw'] * df_output['capacity_mwe'] * 1000 / 1e6
        )
    
    # Use actual construction time if available, otherwise planned
    if 'construction_time_actual' in df_output.columns:
        df_output['construction_time'] = df_output['construction_time_actual'].fillna(
            df_output.get('construction_time_planned', np.nan)
        )
    
    # Remove rows with missing critical data
    critical_cols = ['name', 'capacity_mwe', 'occ_usd_per_kw']
    df_output = df_output.dropna(subset=critical_cols)
    
    # Sort by tech category for easier analysis
    if 'tech_category' in df_output.columns:
        df_output = df_output.sort_values('tech_category')
    
    # Save as semicolon-delimited CSV
    df_output.to_csv(output_csv, sep=';', index=False, encoding='utf-8-sig')
    
    print(f"\n✓ Successfully extracted {len(df_output)} reactors")
    print(f"✓ Saved to: {output_csv}")
    
    # Summary by tech category
    if 'tech_category' in df_output.columns:
        print("\n--- Reactors by Technology Category ---")
        print(df_output.groupby('tech_category').size())
        
        print("\n--- Cost Range by Category (USD/kW) ---")
        cost_summary = df_output.groupby('tech_category')['occ_usd_per_kw'].agg(['count', 'min', 'median', 'max'])
        print(cost_summary)
    
    # Summary by size category
    if 'size_category' in df_output.columns:
        print("\n--- Reactors by Size Category ---")
        print(df_output.groupby('size_category').size())
    
    # Summary by reactor type
    if 'reactor_type' in df_output.columns:
        print("\n--- Reactors by Type ---")
        print(df_output.groupby('reactor_type').size())
    
    print(f"\n--- Columns in Output ---")
    print(list(df_output.columns))
    
    print(f"\n--- Data Completeness ---")
    completeness = df_output.isnull().sum()
    print(completeness[completeness > 0])  # Only show columns with missing data
    
    return df_output


if __name__ == "__main__":
    excel_file = 'reactor_data_raw.xlsx'
    df = extract_reactor_data(excel_file)
    
    if df is not None:
        print("\n--- First Few Rows ---")
        print(df[['name', 'reactor_type', 'size_category', 'tech_category', 'occ_usd_per_kw']].head(10))
        
        # Check for FOAK vs NOAK distribution
        if 'foak_noak' in df.columns and 'tech_category' in df.columns:
            print("\n--- FOAK/NOAK Distribution by Tech Category ---")
            print(pd.crosstab(df['tech_category'], df['foak_noak']))