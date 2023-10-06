import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import shap
from sklearn.inspection import permutation_importance

# # Load SNP data (modify as needed)
# snp_data = np.array([
#     [1.0, 0.0, np.nan, 0.0],
#     [2.0, 0.0, np.nan, 2.0],
#     [0.0, 1.0, 2.0, 0.0]
# ])
#
# # Convert NaN values to a specific value (e.g., -1) to represent missing data
# snp_data[np.isnan(snp_data)] = -1
#
# # Load phenotype labels as a NumPy array (modify as needed)
# phenotype_labels = np.array([
#     [7.],
#     [9.],
#     [2.]
#
# ])

snp_data = np.load('snp2.npy')
snp_data[np.isnan(snp_data)] = -1
phenotype_labels = np.load('pheno1.npy')

#print (snp_data)
#print (phenotype_labels)

# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(snp_data, phenotype_labels, test_size=0.2, random_state=42)

# Standardize the input features (optional but recommended)
scaler = StandardScaler()
X_train = scaler.fit_transform(X_train)
X_test = scaler.transform(X_test)

# Define a deep learning model
model = keras.Sequential([
    layers.Input(shape=(X_train.shape[1],)),
    layers.Dense(128, activation='relu'),
    layers.Dense(64, activation='relu'),
    layers.Dense(1, activation='linear')  # Adjust the output layer for regression or classification
])

# Compile the model
model.compile(optimizer='adam', loss='mean_squared_error')  # Adjust loss function for your task

# Train the model
history = model.fit(X_train, y_train, epochs=50, batch_size=32, validation_split=0.2)

# Evaluate the model on the test set
test_loss = model.evaluate(X_test, y_test)
print(f'Test Loss: {test_loss}')

num_features = X_train.shape[1]

# Ensure the weights are in a 1D array and not a 2D array
first_layer_weights = model.layers[1].get_weights()[0]

first_layer_weights = first_layer_weights.flatten()

# Create a bar chart to visualize the weights
plt.figure(figsize=(10, 6))
plt.bar(range(num_features), first_layer_weights)
plt.xlabel('Input Features')
plt.ylabel('Weight Magnitude')
plt.title('Feature Importance in First Hidden Layer')
plt.xticks(range(num_features), range(1, num_features + 1))  # Adjust for feature labels if needed
plt.show()

# # Calculate permutation feature importances (custom implementation)
# n_permutations = 100  # Number of permutations to perform
# feature_importances = np.zeros(X_test.shape[1])
# print (feature_importances.max())
#
# for _ in range(n_permutations):
#     shuffled_X_test = X_test.copy()
#     np.random.shuffle(shuffled_X_test)  # Permute the feature values
#     shuffled_loss = model.evaluate(shuffled_X_test, y_test, verbose=0)
#     importance = test_loss - shuffled_loss
#     feature_importances += importance
#
# # Normalize feature importances
# feature_importances /= n_permutations
#
# # Print feature importances
# print("Feature Importances:")
# for i, importance in enumerate(feature_importances):
#     print(f"SNP{i+1}: {importance:.4f}")
#
# # Create a Manhattan plot for feature importance
# plt.figure(figsize=(12, 6))
# plt.bar(range(len(feature_importances)), -np.log10(feature_importances), color='royalblue')
# plt.title('Manhattan Plot for Feature Importance')
# plt.xlabel('SNP')
# plt.ylabel('-log10(Importance)')
# plt.xticks(range(len(feature_importances)), [f'SNP{i+1}' for i in range(len(feature_importances))], rotation=90)
# plt.show()
#
# # Make predictions on new data
# new_data = np.array([[0.0, 1.0, 2.0, 0.0]])  # Replace with your new SNP data
# new_data[np.isnan(new_data)] = -1  # Replace NaN values with the same missing value used above
# new_data = scaler.transform(new_data)  # Standardize new data if necessary
# predictions = model.predict(new_data)
# print (predictions)
# # Save the trained model for future use
# model.save('phenotype_prediction_model.h5')
