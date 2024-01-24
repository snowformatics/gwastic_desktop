import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import xgboost as xgb
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.metrics import accuracy_score, mean_squared_error
from fastlmm.inference import FastLMM

class GWASAI:
    def run_random_forest(self, snp_data, pheno_data, snp_ids, test_size, estimators, gwas_result_name, bed_gp, pheno_gp,
                    genomic_predict, genomic_predict_name, model_nr):
        snp_data[np.isnan(snp_data)] = -1

        # Split data into training and testing sets
        X_train, X_test, y_train, y_test = train_test_split(snp_data, pheno_data, test_size=test_size)#, random_state=42)

        # # Standardize the input features
        ### for GP scalling might be not good
        if genomic_predict:
            print ('no scaling')
        else:
            print ('scaling')
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
            for col, score in zip(snp_ids, rf_model.feature_importances_):
                #f.write(str(col) + ' ' + str(score) + '\n')
                data.append((col, score))
            # Convert the list of tuples into a DataFrame
            df = pd.DataFrame(data, columns=['SNP', 'PValue'])
            df[['Chr', 'ChrPos']] = df['SNP'].str.split(':', expand=True)
            return df
            #return sorted_features_df

    def run_xgboost(self, snp_data, pheno_data, snp_ids, test_size, estimators, gwas_result_name, bed_gp, pheno_gp,
                    genomic_predict, genomic_predict_name, max_dep_set, model_nr):

        snp_data[np.isnan(snp_data)] = -1

        # Split data into training and testing sets
        X_train, X_test, y_train, y_test = train_test_split(snp_data, pheno_data, test_size=test_size)#, random_state=42)

        if genomic_predict:
            print ('no scaling')
        else:
            print ('scaling')
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
            for col, score in zip(snp_ids, xgb_model.feature_importances_):
                #f.write(str(col) + ' ' + str(score) + '\n')
                data.append((col, score))
            # Convert the list of tuples into a DataFrame
            df = pd.DataFrame(data, columns=['SNP', 'PValue'])
            df[['Chr', 'ChrPos']] = df['SNP'].str.split(':', expand=True)
            return df

    def run_lmm_gp(self, snp_data, pheno_data, snp_ids, test_size, estimators, gwas_result_name, bed_gp, pheno_gp,
                    genomic_predict, genomic_predict_name, max_dep_set, model_nr):
        import statsmodels.api as sm
        from statsmodels.genmod.bayes_mixed_glm import BinomialBayesMixedGLM
        def calculate_grm(snp_data):
            """ Calculate the Genomic Relationship Matrix (GRM) """
            # Center and scale SNP data
            snp_data = (snp_data - np.mean(snp_data, axis=0)) / np.std(snp_data, axis=0)
            # Calculate GRM
            grm = np.dot(snp_data, snp_data.T) / snp_data.shape[1]
            return grm


        snp_data[np.isnan(snp_data)] = -1

        # Split data into training and testing sets
        X_train, X_test, y_train, y_test = train_test_split(snp_data, pheno_data, test_size=test_size)#, random_state=42)

        # Calculate GRM
        K = calculate_grm(X_train)

        # Add a column of ones to X_train to represent the fixed effects (intercept)
        X_train_with_intercept = sm.add_constant(X_train)

        # Fit the linear mixed model
        model = sm.MixedLM(y_train, X_train_with_intercept, groups=np.ones(len(y_train)), exog_re=np.eye(len(y_train)),
                           re_formula="~0 + x")
        fit = model.fit()
    #
    #     # # Standardize the input features
    #     ### for GP scalling might be not good
    #     if genomic_predict:
    #         print ('no scaling')
    #     else:
    #         print ('scaling')
    #         #scaler = StandardScaler()
    #         #X_train = scaler.fit_transform(X_train)
    #     #X_test = scaler.transform(X_test)
    #     #train = snp_data[:-10, :]
    #     #test = snp_data[-10:, :]
    #
    #     #print (pheno_data)
    #
    #     # Define and train a Random Forest model
    #     fastlmm_model = FastLMM(GB_goal=2)
    #     fastlmm_model.fit(X=snp_data, y=pheno_data)
    #     print ('ok')
    #     #mean, covariance = fastlmm_model.predict(K0_whole_test=test)
    #
    #
    #     "Predicted means and stdevs"
    #     print (mean.val[:, 0])
    #     print (np.sqrt(np.diag(covariance.val)))
    #
    #     # #mean, covariance = fastlmm.predict(K0_whole_test=X_test)
    #     # bed_data = bed_gp.read().iid
    #     # pheno_data = pheno_gp.read().iid
    #     # df_gp = pd.DataFrame(bed_data, columns=['ID1', 'BED_ID2'])
    #     # predicted_values = fastlmm_model.predict(bed_gp.read().val)
    #     # df_gp['Predicted_Value'] = predicted_values
    #     # print (df_gp)
    #
    #

        # # Plot actual phenotype and predicted phenotype
        # whole_pheno = Pheno(pheno_fn)
        # actual_pheno = whole_pheno[whole_pheno.iid_to_index(mean.iid), :].read()
        # pylab.plot(actual_pheno.val, "r.")
        # pylab.plot(mean.val, "b.")
        # pylab.errorbar(np.arange(mean.iid_count), mean.val, yerr=np.sqrt(np.diag(covariance.val)), fmt='.')
        # pylab.xlabel('testing examples')
        # pylab.ylabel('phenotype, actual (red) and predicted (blue with stdev)')
        # pylab.show()



        #
        # if genomic_predict:
        #     bed_data = bed_gp.read().iid
        #     pheno_data = pheno_gp.read().iid
        #     df_gp = pd.DataFrame(bed_data, columns=['ID1', 'BED_ID2'])
        #     predicted_values = rf_model.predict(bed_gp.read().val)
        #     df_gp['Predicted_Value'] = predicted_values
        #     #print (df_gp)
        #     df_pheno = pd.DataFrame(pheno_data, columns=['ID1', 'Pheno_ID2'])
        #     df_pheno['Pheno_Value'] = pheno_gp.read().val
        #
        #     merged_df = df_gp.merge(df_pheno, on='ID1', how='outer')
        #     merged_df['Difference'] = (merged_df['Pheno_Value'] - merged_df['Predicted_Value']).abs()
        #     merged_df['Predicted_Value'] = merged_df['Predicted_Value'].astype(float)
        #     merged_df['Predicted_Value'] = merged_df['Predicted_Value'].round(3)
        #     merged_df['Difference'] = merged_df['Difference'].round(3)
        #     merged_df.to_csv(genomic_predict_name, sep=';', index=False)
        #     return merged_df



