import numpy as np
import matplotlib.pyplot as plt

def question4_pca_reconstruction_loss():
    """
    Question 4: PCA Reconstruction Loss Analysis
    
    This function analyzes class 4 data from the dataset by:
    1. Extracting class 4 samples and mean-centering
    2. Computing PCA using SVD
    3. Calculating reconstruction loss for different numbers of components (L)
    4. Plotting reconstruction loss vs. number of components
    5. Providing additional analysis on variance explained
    """
    
    # Load the dataset
    datafile = np.load('p3378_dr_dataset.npz')
    print("Dataset loaded successfully!")
    
    # --- 1. Extract data for class 4 ---
    print("\n=== Step 1: Extract Class 4 Data ===")
    X_class4 = datafile['meanimgs'][datafile['labels'] == 4]
    N, D = X_class4.shape
    
    print(f"Total samples in dataset: {len(datafile['labels'])}")
    print(f"Number of classes: {len(np.unique(datafile['labels']))}")
    print(f"Class 4 samples: N = {N}")
    print(f"Feature dimensions: D = {D}")
    print(f"Class 4 data shape: {X_class4.shape}")
    
    # --- 2. Mean-center the data ---
    print("\n=== Step 2: Mean-Center the Data ===")
    # PCA requires the data to be centered around the origin
    X_mean = X_class4.mean(axis=0)
    X = X_class4 - X_mean
    
    print(f"Original data mean: {np.mean(X_class4):.6f}")
    print(f"Centered data mean: {np.mean(X):.10f} (should be ~0)")
    print(f"Mean vector shape: {X_mean.shape}")
    
    # --- 3. Compute PCA ---
    print("\n=== Step 3: Compute PCA using SVD ===")
    # SVD is the standard method for finding principal components
    # The rows of Vt are the principal components, sorted by importance
    U, s, Vt = np.linalg.svd(X, full_matrices=False)
    # Transpose Vt to have principal components as columns
    principal_components = Vt.T
    
    print(f"SVD decomposition completed:")
    print(f"  U shape (left singular vectors): {U.shape}")
    print(f"  s shape (singular values): {s.shape}")
    print(f"  Vt shape (right singular vectors): {Vt.shape}")
    print(f"  Principal components shape: {principal_components.shape}")
    
    # Print the first few singular values
    print(f"\nTop 10 singular values:")
    for i in range(min(10, len(s))):
        print(f"  s[{i}] = {s[i]:.4f}")
    
    # --- 4. Loop to compute reconstruction loss ---
    print("\n=== Step 4: Compute Reconstruction Loss ===")
    # This follows the method described in the problem statement
    reconstruction_losses = []
    L_values = range(1, D + 1)
    
    print(f"Computing reconstruction losses for L = 1 to {D}...")
    
    for L in L_values:
        # a. Select the top L principal components
        U_L = principal_components[:, :L]
        
        # b. Reconstruct the data matrix X_hat_L
        # Formula: X_hat_L = (X @ U_L) @ U_L.T
        X_hat_L = (X @ U_L) @ U_L.T
        
        # c. Compute the squared Frobenius norm of the residual
        residual = X - X_hat_L
        loss = np.sum(residual**2)  # Frobenius norm squared
        reconstruction_losses.append(loss)
        
        # Print progress for key values
        if L in [1, 5, 10, 20, 30, 40, 50] or L == D:
            variance_explained = 1 - loss/reconstruction_losses[0] if reconstruction_losses[0] > 0 else 0
            print(f"  L = {L:3d}: Loss = {loss:.2e}, Variance explained = {variance_explained:.4f}")
    
    # Statistical analysis
    print(f"\n=== Reconstruction Loss Statistics ===")
    total_variance = reconstruction_losses[0]  # Loss when L=1 (maximum loss)
    print(f"Maximum loss (L=1): {total_variance:.2e}")
    print(f"Minimum loss (L={D}): {reconstruction_losses[-1]:.2e}")
    print(f"Loss reduction ratio: {reconstruction_losses[-1]/total_variance:.2e}")
    
    # Find L where certain percentages of variance are explained
    print(f"\n=== Variance Explained Analysis ===")
    for threshold in [0.8, 0.9, 0.95, 0.99]:
        target_loss = total_variance * (1 - threshold)
        L_threshold = None
        for i, loss in enumerate(reconstruction_losses):
            if loss <= target_loss:
                L_threshold = i + 1
                break
        if L_threshold:
            print(f"L for {threshold*100:.0f}% variance explained: {L_threshold}")
    
    # --- 5. Plot the results to match the example image ---
    print("\n=== Step 5: Generate Plots ===")
    plt.figure(figsize=(12, 8))
    
    # a. Limit the plot range to match the example's x-axis
    plot_L_max = 50
    if D < plot_L_max:
        plot_L_max = D
    
    # b. Plot the curve with specified markers and style
    plt.plot(L_values[:plot_L_max], reconstruction_losses[:plot_L_max],
             marker='o', markersize=4, linestyle='-', linewidth=2, label='class 4')
    
    # c. Set titles and labels to match the example text
    plt.title('PCA Reconstruction Loss vs. Number of Components', fontsize=14, fontweight='bold')
    plt.xlabel('number of PCA components (L)', fontsize=12)  # Fixed: removed "L2 norm"
    plt.ylabel('reconstruction loss', fontsize=12)
    
    # d. Add the legend
    plt.legend(fontsize=11)
    
    # e. Format the y-axis to use scientific notation with math text
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useMathText=True)
    
    # f. Set x-axis limits to add some margin
    plt.xlim(-2, plot_L_max + 2)
    
    # g. Enable grid for better readability
    plt.grid(True, alpha=0.3)
    
    # h. Improve layout and show the plot
    plt.tight_layout()
    plt.show()
    
    # Additional plot: Explained Variance Ratio
    plt.figure(figsize=(12, 8))
    explained_variance_ratio = [(total_variance - loss) / total_variance for loss in reconstruction_losses]
    
    plt.plot(L_values[:plot_L_max], explained_variance_ratio[:plot_L_max],
             marker='s', markersize=4, linestyle='-', linewidth=2, 
             color='red', label='Explained Variance Ratio')
    
    plt.title('Explained Variance Ratio vs. Number of Components', fontsize=14, fontweight='bold')
    plt.xlabel('number of PCA components (L)', fontsize=12)
    plt.ylabel('explained variance ratio', fontsize=12)
    plt.ylim(0, 1.05)
    plt.xlim(-2, plot_L_max + 2)
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()
    
    return {
        'reconstruction_losses': reconstruction_losses,
        'explained_variance_ratio': explained_variance_ratio,
        'singular_values': s,
        'principal_components': principal_components,
        'N': N,
        'D': D
    }

# Main execution
if __name__ == "__main__":
    print("="*60)
    print("Question 4: PCA Reconstruction Loss Analysis")
    print("="*60)
    
    # Run the analysis
    results = question4_pca_reconstruction_loss()
    
    print(f"\n{'='*60}")
    print("Analysis completed successfully!")
    print(f"Results summary:")
    print(f"- Class 4 samples: {results['N']}")
    print(f"- Feature dimensions: {results['D']}")
    print(f"- Principal components computed: {results['principal_components'].shape[1]}")
    print(f"- Reconstruction losses calculated for L=1 to {results['D']}")
    print(f"{'='*60}")