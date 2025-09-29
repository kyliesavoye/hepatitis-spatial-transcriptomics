import sys
from funcs import *
import squidpy as sq
import os
import seaborn as sns
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu
from itertools import combinations


def main():
    # --- Argument Check ---
    # expect 3 items in total
    #if len(sys.argv) != 3:
     #   print("Error: This script requires exactly two file path arguments.")
      #  print(f"Usage: python {sys.argv[0]} <counts_file> <meta_file>")
       # sys.exit(1) 

    # --- Get Arguments ---
    #counts_file_path = sys.argv[1]
    #meta_file_path = sys.argv[2]

    # --- Print Information ---
    print("--- Python Script Starting ---")
    #print(f"Received Counts File: {counts_file_path}")
    #print(f"Received Metadata File: {meta_file_path}")
    print("------------------------------")

    # --- Dummy Analysis ---
    # This is where your actual data analysis logic would go.
    # For this example, we'll just pretend to do some work.
    try:
        ##########################################
        ##call the definitions here (main logic)##
        ##########################################
        
        
        ### COMPARISON ANALYSIS ###
        
        # --- Get dataset-specific paths for all conditions ---
        output_dirs = {}
        for condition in ["AIH", "D", "SN"]:
            paths = get_paths(condition)
            output_dirs[condition] = paths["outdir"]
        
        # Create comparison output directory
        comparison_paths = get_paths("comparison")
        comparison_outdir = comparison_paths["outdir"]
        
        # # --- Load the average z-score CSV from each condition ---
        # zscore_matrices = {}
        # for condition, folder in output_dirs.items():
        #     csv_path = os.path.join(folder, f"{condition}_neighborhood_enrichment_average_radius50.csv")  
        #     if os.path.exists(csv_path):
        #         zscore_matrices[condition] = pd.read_csv(csv_path, index_col=0)
        #         print(f"[INFO] Loaded {condition} matrix from {csv_path}")
        #     else:
        #         print(f"[WARNING] No CSV found for {condition} in {folder}")

        # # --- Example: Compare D vs AIH ---
        # if "D" in zscore_matrices and "AIH" in zscore_matrices:
        #     diff_df = zscore_matrices["D"] - zscore_matrices["AIH"]
        #     plot_neighbourhood_difference(
        #         diff_df,
        #         title="D vs AIH - Δ Neighborhood Enrichment",
        #         output_path=comparison_outdir + "D_vs_AIH_diff.png"
        #     )

        # # --- Example: Compare D vs SN ---
        # if "SN" in zscore_matrices and "D" in zscore_matrices:
        #     diff_df = zscore_matrices["SN"] - zscore_matrices["D"]
        #     plot_neighbourhood_difference(
        #         diff_df,
        #         title="SN vs D - Δ Neighborhood Enrichment",
        #         output_path=comparison_outdir + "SN_vs_D_diff.png"
        #     )

        # # --- Example: Compare AIH vs SN ---
        # if "SN" in zscore_matrices and "AIH" in zscore_matrices:
        #     diff_df = zscore_matrices["SN"] - zscore_matrices["AIH"]
        #     plot_neighbourhood_difference(
        #         diff_df,
        #         title="SN vs AIH - Δ Neighborhood Enrichment",
        #         output_path=comparison_outdir + "SN_vs_AIH_diff.png"
        #     )
        
        cell_pairs_results = {}
        for condition, folder in output_dirs.items():
            csv_path = os.path.join(folder,f"{condition}_combined_zscore_per_fov.csv")
            if os.path.exists(csv_path):
                cell_pairs_results[condition] = pd.read_csv(csv_path)
                print(f"[INFO] Loaded {condition} cell pairs from {csv_path}")
            else:
                print(f"[WARNING] No cell pairs CSV found for {condition} in {folder}")
        print(cell_pairs_results.keys())
        
        # --- Plot boxplots for specific cell pairs ---
        cell_pairs_of_interest = [
            "Macrophage ~ Macrophage",
            "Hepatocyte ~ Macrophage", 
            "Macrophage ~ T cell",
            "Macrophage ~ Monocyte"
        ]
        
        import matplotlib.pyplot as plt
        
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        axes = axes.flatten()
        
        for i, pair_name in enumerate(cell_pairs_of_interest):
            # Collect data for this cell pair across all conditions
            plot_data = []
            condition_data = {}
         
            for condition in ["AIH", "D", "SN"]:  # Ensure consistent order
                if condition in cell_pairs_results:
                    df = cell_pairs_results[condition]
                    # Filter for the specific cell pair using the 'pair' column
                    mask = df['pair'] == pair_name
                    values = df[mask]['zscore'].dropna()
                    condition_data[condition] = values.tolist()
                    for val in values:
                        plot_data.append({"Condition": condition, "Z-score": val})
                else:
                    condition_data[condition] = []
            
            if plot_data:
                df_plot = pd.DataFrame(plot_data)
                sns.boxplot(data=df_plot, x="Condition", y="Z-score", ax=axes[i], order=["AIH", "D", "SN"])
                axes[i].set_title(f"{pair_name}")
                axes[i].set_ylabel("Z-score")
            
            # Statistical testing
            # First, Kruskal-Wallis test for overall significance
            valid_conditions = [cond for cond, data in condition_data.items() if len(data) > 0]
            if len(valid_conditions) >= 2:
                all_data = [condition_data[cond] for cond in valid_conditions if len(condition_data[cond]) > 0]
            if len(all_data) >= 2:
                try:
                    kruskal_stat, kruskal_p = kruskal(*all_data)
                    
                    # Specific pairwise comparisons using Mann-Whitney U test
                    pairwise_p_values = {}
                    comparisons = [("SN", "D"), ("SN", "AIH"), ("D", "AIH")]
                    
                    for cond1, cond2 in comparisons:
                        if (cond1 in condition_data and cond2 in condition_data and 
                            len(condition_data[cond1]) > 0 and len(condition_data[cond2]) > 0):
                            stat, p_val = mannwhitneyu(condition_data[cond1], condition_data[cond2], alternative='two-sided')
                            pairwise_p_values[f"{cond1} vs {cond2}"] = p_val
                    
                    # # Add p-values as text on the plot
                    # axes[i].text(0.02, 0.98, f"Kruskal-Wallis p = {kruskal_p:.3f}", 
                    # transform=axes[i].transAxes, fontsize=8, verticalalignment='top')
                    
                    # Add pairwise p-values
                    text_lines = [f"{pair}: p = {p:.2e}" for pair, p in pairwise_p_values.items()]
                    axes[i].text(0.02, 0.85, "\n".join(text_lines), 
                    transform=axes[i].transAxes, fontsize=14, verticalalignment='top')
                
                except Exception as e:
                    print(f"[WARNING] Statistical test failed for {pair_name}: {e}")
            
            if plot_data:
                df_plot = pd.DataFrame(plot_data)
                colors = ["#66C2A5", "#FC8D62", "#8DA0CB"]
                sns.boxplot(data=df_plot, x="Condition", y="Z-score", ax=axes[i], order=["AIH", "D", "SN"], palette=colors)
                axes[i].set_title(f"{pair_name}", fontsize=16)
                axes[i].set_ylabel("Z-score", fontsize=14)
                axes[i].set_xlabel("Condition", fontsize=14)
                axes[i].tick_params(axis='both', which='major', labelsize=12)
            else:
                axes[i].text(0.5, 0.5, f"No data for {pair_name}", 
                    transform=axes[i].transAxes, ha='center', va='center', fontsize=14)
                axes[i].set_title(f"{pair_name}", fontsize=16)
        
        plt.tight_layout()
        plt.savefig(comparison_outdir + "macrophage_cell_pairs_boxplots_with_stats.png", dpi=300, bbox_inches='tight')
        plt.show()
        print("[INFO] Boxplots with statistical tests saved to macrophage_cell_pairs_boxplots_with_stats.png")


    except FileNotFoundError as e:
        print(f"Error: A file was not found.")
        print(e)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred during analysis: {e}")
        sys.exit(1)

    print("--- Python Script Finished Successfully ---")

if __name__ == "__main__":
    main()
