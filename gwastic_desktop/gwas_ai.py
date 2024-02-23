import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import xgboost as xgb
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.metrics import accuracy_score, mean_squared_error
from fastlmm.inference import FastLMM


class GWASAI:
    def run_random_forest(self, snp_data, pheno_data, df_bim, test_size, estimators, gwas_result_name, bed_gp, pheno_gp,
                    genomic_predict, genomic_predict_name, model_nr):
        snp_data[np.isnan(snp_data)] = -1

        # Split data into training and testing sets
        X_train, X_test, y_train, y_test = train_test_split(snp_data, pheno_data, test_size=test_size)#, random_state=42)

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
        rf_model = RandomForestRegressor(n_estimators=estimators, n_jobs=-1)#, random_state=0) #42
        #rf_model = RandomForestClassifier(n_estimators=estimators,max_depth=2, random_state=0)

        # Fit the model to the training data
        rf_model.fit(X_train, y_train)

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

    # def run_lmm_gp(self, snp_data, pheno_data, gwas_result_name, bed_gp, pheno_gp,
    #                  model_nr):
    #     from sklearn.model_selection import KFold
    #     kf = KFold(n_splits=model_nr + 1)
    #     print (pheno_data)
    #     results = []
    #     #counter = 1
    #     # Start cross-validation
    #     for train_indices, test_indices in kf.split(range(snp_data.iid_count)):
    #         # Split data into training and testing
    #
    #         train = snp_data[train_indices, :]
    #         test = snp_data[test_indices, :]
    #
    #         # Train the model
    #         fastlmm = FastLMM(GB_goal=2)
    #         fastlmm.fit(K0_train=train, y=pheno_data)
    #
    #         # Test the model
    #         mean, covariance = fastlmm.predict(K0_whole_test=snp_data)
    #         df = pd.DataFrame({'Column1': mean.val[:, 0], 'Column2': mean.iid[:, 0]})
    #         #df.to_csv(str(counter) + '_lmm_gp_priming2.csv', index=False)
    #         #counter += 1
    #         print(df)
    #
    #         bed_data = bed_gp.read().iid
    #         pheno_data2 = pheno_gp.read().iid
    #         df_gp = pd.DataFrame(bed_data, columns=['ID1', 'BED_ID2'])
    #         predicted_values = df['Column1']
    #         df_gp['Predicted_Value'] = predicted_values
    #
    #         df_pheno = pd.DataFrame(pheno_data2, columns=['ID1', 'Pheno_ID2'])
    #         df_pheno['Pheno_Value'] = pheno_gp.read().val
    #         merged_df = df_gp.merge(df_pheno, on='ID1', how='outer')
    #         merged_df['Difference'] = (merged_df['Pheno_Value'] - merged_df['Predicted_Value']).abs()
    #         merged_df['Predicted_Value'] = merged_df['Predicted_Value'].astype(float)
    #         merged_df['Predicted_Value'] = merged_df['Predicted_Value'].round(5)
    #         merged_df['Difference'] = merged_df['Difference'].round(5)
    #         results.append(merged_df)
    #         print (merged_df)
    #     return results
    #
