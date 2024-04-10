import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import xgboost as xgb
from sklearn.ensemble import RandomForestRegressor

###This file is depreciated
class GWASAI:
    """This file is depreciated."""
    def run_random_forest(self, snp_data, pheno_data, df_bim, test_size, estimators, gwas_result_name, bed_gp, pheno_gp,
                    genomic_predict, genomic_predict_name, model_nr):
        print(np.sum(np.isnan(snp_data)))
        snp_data[np.isnan(snp_data)] = -1
        print(np.sum(np.isnan(snp_data)))
        print(snp_data.shape, pheno_data.shape)
        # Split data into training and testing sets
        X_train, X_test, y_train, y_test = train_test_split(snp_data, pheno_data, test_size=test_size, random_state=42)#, random_state=42)
        print(X_train.shape, X_test.shape, y_train.shape, y_test.shape)
        # # Standardize the input features
        ### for GP scalling might be not good
        if genomic_predict:
            pass
        else:
            #print ('scaling')
            scaler = StandardScaler()
            X_train = scaler.fit_transform(X_train)
        #X_test = scaler.transform(X_test)

        # Define and train a Random Forest model
        rf_model = RandomForestRegressor(n_estimators=estimators, n_jobs=-1, random_state=42)#, random_state=0) #42
        #rf_model = RandomForestClassifier(n_estimators=estimators,max_depth=2, random_state=0)
        #print(X_train, y_train, y_test)
        #print (type(X_train), type(y_train))
        #import matplotlib.pyplot as plt
        #import seaborn as sns
        #plt.imshow(X_train, cmap='viridis')  # 'cmap' denotes the color map
        #plt.colorbar()  # Optional: adds a color bar to the side
        #plt.show()
        #df = pd.DataFrame(X_train)
        # Create a heatmap
        #sns.heatmap(df)
        #plt.show()

        # Fit the model to the training data
        rf_model.fit(X_train, y_train)
        print(rf_model.feature_importances_)
        # import seaborn as sns
        # print(len(snp_data))
        #
        # importances = rf_model.feature_importances_
        # indices = np.argsort(importances)[::-1]
        #
        # # Get the indices of the top 10 features
        # top_10_indices = np.argsort(importances)[-50:]
        #
        # # Create a DataFrame with the top 10 SNPs
        # top_10_snps = pd.DataFrame(snp_data[:, top_10_indices], columns=[f'SNP_{i}' for i in top_10_indices])
        #
        # # Calculate the correlation matrix
        # corr = top_10_snps.corr()
        #
        # # Generate a heatmap
        # plt.figure(figsize=(10, 8))
        # sns.heatmap(corr, annot=True, cmap='coolwarm', fmt=".2f", xticklabels=1, yticklabels=1)
        # plt.title("Correlation Heatmap of Top 10 SNPs")
        # plt.show()

        if genomic_predict:
            bed_data = bed_gp.read().iid
            pheno_data = pheno_gp.read().iid
            df_gp = pd.DataFrame(bed_data, columns=['ID1', 'BED_ID2'])
            predicted_values = rf_model.predict(bed_gp.read().val)
            df_gp['Predicted_Value'] = predicted_values
            #print (df_gp)
            df_pheno = pd.DataFrame(pheno_data, columns=['ID1', 'Pheno_ID2'])
            df_pheno['Pheno_Value'] = pheno_gp.read().val

            merged_df = df_gp.merge(df_pheno, on='ID1', how='outer')
            merged_df['Difference'] = (merged_df['Pheno_Value'] - merged_df['Predicted_Value']).abs()
            merged_df['Predicted_Value'] = merged_df['Predicted_Value'].astype(float)
            merged_df['Predicted_Value'] = merged_df['Predicted_Value'].round(3)
            merged_df['Difference'] = merged_df['Difference'].round(3)
            merged_df.to_csv(genomic_predict_name, sep=';', index=False)
            return merged_df

        else:

            # features_dict = dict(zip(snp_ids, rf_model.feature_importances_))
            # # Sort the dictionary by importance
            # sorted_features = sorted(features_dict.items(), key=lambda x: x[1], reverse=True)
            # # Convert to a DataFrame for easy saving
            # sorted_features_df = pd.DataFrame(sorted_features, columns=['Feature', 'Importance'])
            # # Save to a CSV file
            # sorted_features_df.to_csv(gwas_result_name, index=False)
            # sorted_features_df.columns = ['SNP', 'PValue']
            # sorted_features_df[['Chr', 'ChrPos']] = sorted_features_df['SNP'].str.split(':', expand=True)

            # f = open(gwas_result_name, 'w')

            data = []
            snp_ids = df_bim.iloc[:, 1].tolist()
            for col, score in zip(snp_ids, rf_model.feature_importances_):
                data.append((col, score))
            # Convert the list of tuples into a DataFrame
            df = pd.DataFrame(data, columns=['SNP', 'PValue'])
            df = pd.merge(df, df_bim, on='SNP', how='left')
            del df['NA']
            #df[['Chr', 'ChrPos']] = df['SNP'].str.split(':', expand=True)
            return df

    def run_xgboost(self, snp_data, pheno_data, df_bim, test_size, estimators, gwas_result_name, bed_gp, pheno_gp,
                    genomic_predict, genomic_predict_name, max_dep_set, model_nr):

        snp_data[np.isnan(snp_data)] = -1

        # Split data into training and testing sets
        X_train, X_test, y_train, y_test = train_test_split(snp_data, pheno_data, test_size=test_size)#, random_state=42)

        if genomic_predict:
            pass
        else:

            scaler = StandardScaler()
            X_train = scaler.fit_transform(X_train)
        # # Standardize the input features
        #scaler = StandardScaler()
        #X_train = scaler.fit_transform(X_train)
        #X_test = scaler.transform(X_test)

        # Define and train a Random Forest model
        xgb_model = xgb.XGBRegressor(n_estimators=estimators, learning_rate=0.1, max_depth=max_dep_set)#, random_state=0)
        #xgb_model = xgb.XGBRegressor(n_estimators=estimators, random_state=42)

        # Fit the model to the training data
        xgb_model.fit(X_train, y_train)

        # test_predictions = rf_model.predict(X_test)
        # test_loss = np.mean((test_predictions - y_test)**2)  # Mean squared error
        # print(f'Test Loss: {test_loss}')

        if genomic_predict:
            bed_data = bed_gp.read().iid
            pheno_data = pheno_gp.read().iid
            df_gp = pd.DataFrame(bed_data, columns=['ID1', 'BED_ID2'])
            predicted_values = xgb_model.predict(bed_gp.read().val)
            df_gp['Predicted_Value'] = predicted_values

            df_pheno = pd.DataFrame(pheno_data, columns=['ID1', 'Pheno_ID2'])
            df_pheno['Pheno_Value'] = pheno_gp.read().val
            merged_df = df_gp.merge(df_pheno, on='ID1', how='outer')
            merged_df['Difference'] = (merged_df['Pheno_Value'] - merged_df['Predicted_Value']).abs()
            merged_df['Predicted_Value'] = merged_df['Predicted_Value'].astype(float)
            merged_df['Predicted_Value'] = merged_df['Predicted_Value'].round(3)
            merged_df['Difference'] = merged_df['Difference'].round(3)

            merged_df.to_csv(genomic_predict_name, sep=';', index=False)
            return merged_df

        else:
            # features_dict = dict(zip(snp_ids, xgb_model.feature_importances_))
            # # Sort the dictionary by importance
            # sorted_features = sorted(features_dict.items(), key=lambda x: x[1], reverse=True)
            # # Convert to a DataFrame for easy saving
            # sorted_features_df = pd.DataFrame(sorted_features, columns=['Feature', 'Importance'])
            # # Save to a CSV file
            # sorted_features_df.to_csv(gwas_result_name, index=False)
            # sorted_features_df.columns = ['SNP', 'PValue']
            # sorted_features_df[['Chr', 'ChrPos']] = sorted_features_df['SNP'].str.split(':', expand=True)

            #f = open(gwas_result_name, 'w')
            data = []
            snp_ids = df_bim.iloc[:, 1].tolist()
            for col, score in zip(snp_ids, xgb_model.feature_importances_):
                #f.write(str(col) + ' ' + str(score) + '\n')
                data.append((col, score))
            # Convert the list of tuples into a DataFrame
            df = pd.DataFrame(data, columns=['SNP', 'PValue'])
            df = pd.merge(df, df_bim, on='SNP', how='left')
            del df['NA']
            #df[['Chr', 'ChrPos']] = df['SNP'].str.split(':', expand=True)
            return df

