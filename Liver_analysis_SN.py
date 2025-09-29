import sys
from funcs import *
import squidpy as sq
import matplotlib.colors as mcolors


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


        ### ANALYSING SN SAMPLES ###

        # --- Get dataset-specific paths ---
        paths = get_paths("SN")
        counts = paths["counts"]
        meta = paths["meta"]
        outdirectory = paths["outdir"]

        # Load SN samples
        adata_SN_606 = sq.read.nanostring(
            path=".",  # current directory
            counts_file=counts + "SN_606_expression_mat_fixed.csv", 
            meta_file=meta + "SN_606_spatial_metadata_fixed.csv")

        adata_SN_606.obs['cell_type'] = adata_SN_606.obs['final_cell_types'].astype('category')

        adata_SN_1042 = sq.read.nanostring(
            path=".",  # current directory
            counts_file=counts + "SN_1042_expression_mat_fixed.csv",  
            meta_file=meta + "SN_1042_spatial_metadata_fixed.csv"
        )
        adata_SN_1042.obs['cell_type'] = adata_SN_1042.obs['final_cell_types'].astype('category')

        adata_SN_2739 = sq.read.nanostring(
            path=".",  # current directory
            counts_file=counts + "SN_2739_expression_mat_fixed.csv",
            meta_file=meta + "SN_2739_spatial_metadata_fixed.csv"
        )
        adata_SN_2739.obs['cell_type'] = adata_SN_2739.obs['final_cell_types'].astype('category')

        ### Run neighbourhood enrichment for each SN sample at different radii ###

        z606_df_30 = neighbourhood_enrichment(adata_SN_606, "SN_606", radius=30, save_plots=True, outdir=outdirectory)
        z1042_df_30 = neighbourhood_enrichment(adata_SN_1042, "SN_1042", radius=30, save_plots=True, outdir=outdirectory)
        z2739_df_30 = neighbourhood_enrichment(adata_SN_2739, "SN_2739", radius=30, save_plots=True, outdir=outdirectory)

        z_avg_SN_30 = average_neighbourhood_plot(
            [z606_df_30, z1042_df_30, z2739_df_30],
            output_path_png=outdirectory + "SN_neighborhood_enrichment_average_radius30.png",
            output_path_csv=outdirectory + "SN_neighborhood_enrichment_average_radius30.csv",
            title="Average Neighborhood Enrichment SN Samples - Radius 30 microns"
        )

        z606_df_80 = neighbourhood_enrichment(adata_SN_606, "SN_606", radius=80, save_plots=True, outdir=outdirectory)
        z1042_df_80 = neighbourhood_enrichment(adata_SN_1042, "SN_1042", radius=80, save_plots=True, outdir=outdirectory)
        z2739_df_80 = neighbourhood_enrichment(adata_SN_2739, "SN_2739", radius=80, save_plots=True, outdir=outdirectory)

        z_avg_SN_80 = average_neighbourhood_plot(
            [z606_df_80, z1042_df_80, z2739_df_80],
            output_path_png=outdirectory + "SN_neighborhood_enrichment_average_radius80.png",
            output_path_csv=outdirectory + "SN_neighborhood_enrichment_average_radius80.csv",
            title="Average Neighborhood Enrichment SN Samples - Radius 80 microns"
        )

        # Run neighbourhood enrichment for radius 50 microns last as the current 
        # neighbourhood enrichment function modifies the adata object when run

        z606_df_50 = neighbourhood_enrichment(adata_SN_606, "SN_606", radius=50, save_plots=True, outdir=outdirectory)
        z1042_df_50 = neighbourhood_enrichment(adata_SN_1042, "SN_1042", radius=50, save_plots=True, outdir=outdirectory)
        z2739_df_50 = neighbourhood_enrichment(adata_SN_2739, "SN_2739", radius=50, save_plots=True, outdir=outdirectory)

        z_avg_SN_50 = average_neighbourhood_plot(
            [z606_df_50, z1042_df_50, z2739_df_50],
            output_path_png=outdirectory + "SN_neighborhood_enrichment_average_radius50.png",
            output_path_csv=outdirectory + "SN_neighborhood_enrichment_average_radius50.csv",
            title="Average Neighborhood Enrichment - SN Samples"
        )

        ### Set colourmap ###

        # Optional: shared consistent palette across all samples
        cell_types = ["B cell", "Epithelial", "Hepatic stellate", "Hepatocyte", "Macrophage", "Monocyte", "NK cell", "Neutrophil", "Plasma cell", "T cell", "Type 1 LSEC", "Type 2 LSEC"]
        
        # Custom hex palette
        palette = [
            "#FFA6C9",  # B cell
            "#D5E8D4",  # Epithelial
            "#FFD1DC",  # Hepatic stellate
            "#AEC6CF",  # Hepatocyte
            "#F7B2AB",  # Macrophage
            "#FDFD96",  # Monocyte
            "#FFDAC1",  # NK cell
            "#99EE99",  # Neutrophil
            "#A9E2F3",  # Plasma cell
            "#B39EB5",  # T cell
            "#F4BBFF",  # Type 1 LSEC
            "#E6E6FA"   # Type 2 LSEC
        ]

        celltype_color_map = dict(zip(cell_types, palette))

        ### Run Ripley's G functions for each SN sample ###

        sq.gr.ripley(adata_SN_606, cluster_key='cell_type', mode='G', max_dist=100)
        # Make sure cell_type is categorical
        adata_SN_606.obs["cell_type"] = adata_SN_606.obs["cell_type"].astype("category")

        categories = adata_SN_606.obs['cell_type'].cat.categories
        palette_ordered = [celltype_color_map[ct] for ct in categories]
        cmap = mcolors.ListedColormap(palette_ordered)

        sq.pl.ripley(
            adata_SN_606,
            cluster_key='cell_type',
            mode='G',
            palette=cmap,
            plot_sims=True,
            save=outdirectory + "SN_606_ripley_G_sq.png"
        )
        
        sq.gr.ripley(adata_SN_1042, cluster_key='cell_type', mode='G', max_dist=100)

        # Make sure cell_type is categorical
        adata_SN_1042.obs["cell_type"] = adata_SN_1042.obs["cell_type"].astype("category")

        categories = adata_SN_1042.obs['cell_type'].cat.categories
        palette_ordered = [celltype_color_map[ct] for ct in categories]
        cmap = mcolors.ListedColormap(palette_ordered)

        sq.pl.ripley(
            adata_SN_1042,
            cluster_key='cell_type',
            mode='G',
            palette=cmap,
            plot_sims=True,
            save=outdirectory + "SN_1042_ripley_G_sq.png"
        )        
        
        sq.gr.ripley(adata_SN_2739, cluster_key='cell_type', mode='G', max_dist=100)
        # Make sure cell_type is categorical
        adata_SN_2739.obs["cell_type"] = adata_SN_2739.obs["cell_type"].astype("category")

        categories = adata_SN_2739.obs['cell_type'].cat.categories
        palette_ordered = [celltype_color_map[ct] for ct in categories]
        cmap = mcolors.ListedColormap(palette_ordered)

        sq.pl.ripley(
            adata_SN_2739,
            cluster_key='cell_type',
            mode='G',
            palette=cmap,
            plot_sims=True,
            save=outdirectory + "SN_2739_ripley_G_sq.png"
        )  

        # Computing and plotting the average G statistic
        ripleyG_606 = adata_SN_606.uns["cell_type_ripley_G"]
        ripleyG_2739 = adata_SN_2739.uns["cell_type_ripley_G"]
        ripleyG_1042 = adata_SN_1042.uns["cell_type_ripley_G"]

        sims_606G = process_sims_stat(ripleyG_606, '606')
        sims_2739G = process_sims_stat(ripleyG_2739, '2739')
        sims_1042G = process_sims_stat(ripleyG_1042, '1042')

        sims_allG = pd.concat([sims_606G, sims_2739G, sims_1042G], axis=0)

        plt.figure(figsize=(12, 8))

        df_606G = ripleyG_606['G_stat']
        df_2739G = ripleyG_2739['G_stat']
        df_1042G = ripleyG_1042['G_stat']

        df_606G['source'] = '606'
        df_2739G['source'] = '2739'
        df_1042G['source'] = '1042'

        df_allG = pd.concat([df_606G, df_2739G, df_1042G], axis=0)

        df_summaryG = df_allG.groupby(['bins', 'cell_type']).agg(
            mean_stats=('stats', 'mean'),
            lower_stats=('stats', lambda x: np.percentile(x, 5)),
            upper_stats=('stats', lambda x: np.percentile(x, 95))
        ).reset_index()

        sims_summaryG = sims_allG.groupby(['bins']).agg(
            lower_sims=('stats', lambda x: np.percentile(x, 5)),
            upper_sims=('stats', lambda x: np.percentile(x, 95)),
            mean_sims=('stats', 'mean')
        ).reset_index()
        
        for cell_type in df_summaryG['cell_type'].unique():
            subset = df_summaryG[df_summaryG['cell_type'] == cell_type]
            color = celltype_color_map.get(cell_type, (0.5, 0.5, 0.5))
            
            plt.plot(subset['bins'], subset['mean_stats'], label=cell_type, color=color)
            plt.fill_between(subset['bins'], 
                     np.clip(subset['lower_stats'], 0, 1), 
                     np.clip(subset['upper_stats'], 0, 1), 
                     color=color, alpha=0.3)
        
        # Plot the average random expectation line (same for all cell types)
        plt.plot(sims_summaryG['bins'], sims_summaryG['mean_sims'], color='black', linestyle='--', alpha=0.5, label='Random expectation')
        plt.xlabel('Distance bins', fontsize=16, weight='bold')
        plt.ylabel('Average Ripley G statistic', fontsize=16, weight='bold')
        plt.title('Average Ripley G function with Random Expectation', fontsize=18, weight='bold')
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12)
        plt.tight_layout()
        plt.savefig(outdirectory + "SN_average_ripley_G_function_with_random_expectation.png")
        plt.close()
        
        # now for selected cell types (Immune cells)
        # Define the subset of cell types you want to plot
        selected_cell_types = ['B cell', 'T cell', 'Macrophage', 'Monocyte', 'NK cell', 'Neutrophil', "Plasma cell"]
        
        plt.figure(figsize=(12, 8))
        
        for cell_type in selected_cell_types:
            subset = df_summaryG[df_summaryG['cell_type'] == cell_type]
            color = celltype_color_map.get(cell_type, (0.5, 0.5, 0.5))
            
            plt.plot(subset['bins'], subset['mean_stats'], label=cell_type, color=color)
            plt.fill_between(subset['bins'], 
                     np.clip(subset['lower_stats'], 0, 1), 
                     np.clip(subset['upper_stats'], 0, 1), 
                     color=color, alpha=0.3)
        
        # Plot the average random expectation line (same for all cell types)
        plt.plot(sims_summaryG['bins'], sims_summaryG['mean_sims'], color='black', linestyle='--', alpha=0.5, label='Random expectation')

        plt.xlabel('Distance bins', fontsize=16, weight='bold')
        plt.ylabel('Average Ripley G statistic', fontsize=16, weight='bold')
        plt.title('Average Ripley G function with Random Expectation (Immune cells)', fontsize=18, weight='bold')
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12)
        plt.tight_layout()
        plt.savefig(outdirectory + "SN_average_ripley_G_immunecells.png")
        plt.close()
        
        # Computing P vs NP G statistic
        adata_SN_1042.obs['P_VS_NP'] = adata_SN_1042.obs["P_VS_NP"].astype("category")
        sq.gr.ripley(adata_SN_1042,cluster_key="P_VS_NP",mode="G", n_neigh=5, n_simulations=100, max_dist=100, n_steps=50, copy=False)
        sq.pl.ripley(adata_SN_1042, cluster_key="P_VS_NP", mode ='G', plot_sims=True, save = outdirectory + "SN_1042_ripley_G_sq_PvsNP.png")
        
        adata_SN_606.obs['P_VS_NP'] = adata_SN_606.obs["P_VS_NP"].astype("category")
        sq.gr.ripley(adata_SN_606,cluster_key="P_VS_NP",mode="G", n_neigh=5, n_simulations=100, max_dist=100, n_steps=50, copy=False)
        sq.pl.ripley(adata_SN_606, cluster_key="P_VS_NP", mode ='G', plot_sims=True, save = outdirectory + "SN_606_ripley_G_sq_PvsNP.png")
        
        adata_SN_2739.obs['P_VS_NP'] = adata_SN_2739.obs["P_VS_NP"].astype("category")
        sq.gr.ripley(adata_SN_2739,cluster_key="P_VS_NP",mode="G", n_neigh=5, n_simulations=100, max_dist=100, n_steps=50, copy=False)
        sq.pl.ripley(adata_SN_2739, cluster_key="P_VS_NP", mode ='G', plot_sims=True, save = outdirectory + "SN_2739_ripley_G_sq_PvsNP.png")
        
        cell_types_P_NP = ['P', 'NP']
        
        palette_2= [(0.9678, 0.4413, 0.5358),(0.9036, 0.5120, 0.1959)]
        
        celltype_color_map_P_NP = dict(zip(cell_types_P_NP, palette_2))
        
        ripleyG_606_PNP = adata_SN_606.uns["P_VS_NP_ripley_G"]
        ripleyG_2739_PNP = adata_SN_2739.uns["P_VS_NP_ripley_G"]
        ripleyG_1042_PNP = adata_SN_1042.uns["P_VS_NP_ripley_G"]

        plt.figure(figsize=(12, 8))

        df_606G = ripleyG_606_PNP['G_stat']
        df_2739G = ripleyG_2739_PNP['G_stat']
        df_1042G = ripleyG_1042_PNP['G_stat']

        df_606G['source'] = '606'
        df_2739G['source'] = '2739'
        df_1042G['source'] = '1042'


        df_allG = pd.concat([df_606G, df_2739G, df_1042G], axis=0)

        df_summaryG = df_allG.groupby(['bins', 'P_VS_NP']).agg(
            mean_stats=('stats', 'mean'),
            lower_stats=('stats', lambda x: np.percentile(x, 5)),
            upper_stats=('stats', lambda x: np.percentile(x, 95))
        ).reset_index()

        sims_summaryG = sims_allG.groupby(['bins']).agg(
            lower_sims=('stats', lambda x: np.percentile(x, 5)),
            upper_sims=('stats', lambda x: np.percentile(x, 95)),
            mean_sims=('stats', 'mean')
        ).reset_index()
        
        for cell_type in df_summaryG['P_VS_NP'].unique():
            subset = df_summaryG[df_summaryG['P_VS_NP'] == cell_type]
            color = celltype_color_map_P_NP.get(cell_type, (0.5, 0.5, 0.5))
            
            plt.plot(subset['bins'], subset['mean_stats'], label=cell_type, color=color)
            plt.fill_between(subset['bins'], 
                     np.clip(subset['lower_stats'], 0, 1), 
                     np.clip(subset['upper_stats'], 0, 1), 
                     color=color, alpha=0.3)
        
        # Plot the average random expectation line (same for all cell types)
        plt.plot(sims_summaryG['bins'], sims_summaryG['mean_sims'], color='black', linestyle='--', alpha=0.5, label='Random expectation')
        plt.xlabel('Distance bins', fontsize=16, weight='bold')
        plt.ylabel('Average Ripley G statistic', fontsize=16, weight='bold')
        plt.title('Average Ripley G function with Random Expectation (P vs NP)', fontsize=18, weight='bold')
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=16)
        plt.tight_layout()
        plt.savefig(outdirectory + "SN_average_ripley_G_function_with_random_expectation_P_VS_NP.png")
        plt.close()

        ### Now run neighbourhood enrichment per fov ###

        z_606_per_fov = run_neighbourhood_enrichment_per_fov(adata_SN_606, "SN_606")
        z_1042_per_fov = run_neighbourhood_enrichment_per_fov(adata_SN_1042, "SN_1042")
        z_2739_per_fov = run_neighbourhood_enrichment_per_fov(adata_SN_2739, "SN_2739")

        # Plot boxplots for each FOV
        z_606_per_fov['pair'] = z_606_per_fov.apply(canonical_pair, axis=1)

        # Compute median z-score for sorting
        pair_order = (
            z_606_per_fov.groupby('pair')['zscore']
            .median()
            .sort_values(ascending=False)
            .index.tolist()
        )
        
        top20_pairs=pair_order[:20]
        # Optional: drop the original columns if not needed
        # z_all = z_all.drop(columns=['cell_type_1', 'cell_type_2'])

        plt.figure(figsize=(20, 8))

        sns.boxplot(data=z_606_per_fov[z_606_per_fov['pair'].isin(top20_pairs)], 
            x="pair", 
            y="zscore", 
            hue="sample", 
            palette=["#66C2A5"], 
            order=top20_pairs
        )

        plt.xticks(rotation=90, fontsize=18)   
        plt.yticks(fontsize=18)              
        plt.xlabel("Cell-Type Pairs", fontsize=18, weight="bold")
        plt.ylabel("Z-Score", fontsize=18, weight="bold")
        plt.title("Top 20 Neighborhood Enrichment Z-Scores per Cell-Type Pair in SN 606 (Ranked by Median)",fontsize=20, weight="bold")
        plt.xlabel("Cell-Type Pairs")
        plt.ylabel("Z-Score")
        plt.legend(title="Sample", bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=15, title_fontsize=16)
        plt.tight_layout()
        plt.savefig(outdirectory + "top20_SN_606_neighborhood_enrichment_zscores_per_pair.png", dpi=300)
        print('top 20 figure saved')
        plt.close()

        z_1042_per_fov['pair'] = z_1042_per_fov.apply(canonical_pair, axis=1)

        pair_order_1042 = (
            z_1042_per_fov.groupby('pair')['zscore']
            .median()
            .sort_values(ascending=False)
            .index.tolist()
        )
        top20_pairs=pair_order_1042[:20]

        plt.figure(figsize=(20, 8))
        sns.boxplot(data=z_1042_per_fov[z_1042_per_fov['pair'].isin(top20_pairs)], 
            x="pair", 
            y="zscore", 
            hue="sample", 
            palette=["#FC8D62"], 
            order=top20_pairs
        )

        plt.xticks(rotation=90, fontsize=18)  
        plt.yticks(fontsize=18)               
        plt.xlabel("Cell-Type Pairs", fontsize=18, weight="bold")
        plt.ylabel("Z-Score", fontsize=18, weight="bold")
        plt.title("Top 20 Neighborhood Enrichment Z-Scores per Cell-Type Pair in SN 1042 (Ranked by Median)",fontsize=20, weight="bold")
        plt.xlabel("Cell-Type Pairs")
        plt.ylabel("Z-Score")
        plt.legend(title="Sample", bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=15, title_fontsize=16)
        plt.tight_layout()
        plt.savefig(outdirectory + "top20_SN_1042_neighborhood_enrichment_zscores_per_pair.png", dpi=300)
        print('top 20 figure saved')
        plt.close()


        z_2739_per_fov['pair'] = z_2739_per_fov.apply(canonical_pair, axis=1)

        pair_order_2739 = (
            z_2739_per_fov.groupby('pair')['zscore']
            .median()
            .sort_values(ascending=False)
            .index.tolist()
        )
        top20_pairs=pair_order_2739[:20]

        plt.figure(figsize=(20, 8))
        sns.boxplot(data=z_2739_per_fov[z_2739_per_fov['pair'].isin(top20_pairs)], 
            x="pair", 
            y="zscore", 
            hue="sample", 
            palette=["#8DA0CB"], 
            order=top20_pairs
        )

        plt.xticks(rotation=90, fontsize=18)   
        plt.yticks(fontsize=18)                
        plt.xlabel("Cell-Type Pairs", fontsize=18, weight="bold")
        plt.ylabel("Z-Score", fontsize=18, weight="bold")
        plt.title("Top 20 Neighborhood Enrichment Z-Scores per Cell-Type Pair in SN 2739 (Ranked by Median)",fontsize=20, weight="bold")
        plt.xlabel("Cell-Type Pairs")
        plt.ylabel("Z-Score")
        plt.legend(title="Sample", bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=15, title_fontsize=16)
        plt.tight_layout()
        plt.savefig(outdirectory + "top20_SN_2739_neighborhood_enrichment_zscores_per_pair.png", dpi=300)
        print('top 20 figure saved')
        plt.close()

        
        # Combine the datasets
        combined_data = pd.concat([z_2739_per_fov, z_606_per_fov, z_1042_per_fov], ignore_index=True)

        # Compute median z-score across both samples for consistent ordering
        pair_order = (
            combined_data.groupby('pair')['zscore']
            .median()
            .sort_values(ascending=False)
            .index.tolist()
        )
        top20_pairs=pair_order[:20]

        # Create the joint boxplot
        plt.figure(figsize=(24, 10))

        sns.boxplot(
            data=combined_data[combined_data['pair'].isin(top20_pairs)], 
            x="pair", 
            y="zscore", 
            hue="sample", 
            palette=["#8DA0CB","#66C2A5", "#FC8D62"], 
            order=top20_pairs
        )

        plt.xticks(rotation=90, fontsize=18)  
        plt.yticks(fontsize=18)                
        plt.xlabel("Cell-Type Pairs", fontsize=18, weight="bold")
        plt.ylabel("Z-Score", fontsize=18, weight="bold")
        plt.title("Top 20 Neighborhood Enrichment Z-Scores per Cell-Type Pair: SN Samples (Ranked by Median)",fontsize=20, weight="bold")
        plt.xlabel("Cell-Type Pairs")
        plt.ylabel("Z-Score")
        plt.legend(title="Sample", bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=15, title_fontsize=16)
        plt.tight_layout()
        plt.savefig(outdirectory + "top20_SN_joint_neighborhood_enrichment_zscores_per_pair.png", dpi=300, bbox_inches='tight')
        print('top 20 figure saved')
        plt.close()
       
        # Combine the datasets
        combined_data = pd.concat([z_606_per_fov, z_1042_per_fov, z_2739_per_fov], ignore_index=True)

        # Save combined z-scores as CSV
        combined_data.to_csv(outdirectory + "SN_combined_zscore_per_fov.csv", index=False)
          
        
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
