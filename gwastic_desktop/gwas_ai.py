import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import xgboost as xgb
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.metrics import accuracy_score, mean_squared_error


class GWASAI:
    def run_random_forest(self, snp_data, pheno_data, snp_ids, test_size, estimators, gwas_result_name, predict):
        print(predict)
        snp_data[np.isnan(snp_data)] = -1

        # Split data into training and testing sets
        X_train, X_test, y_train, y_test = train_test_split(snp_data, pheno_data, test_size=test_size, random_state=42)

        # # Standardize the input features
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train)
        X_test = scaler.transform(X_test)

        # Define and train a Random Forest model
        rf_model = RandomForestRegressor(n_estimators=estimators, random_state=42)

        # Fit the model to the training data
        rf_model.fit(X_train, y_train)

        # Evaluate the model on the test set
        #test_predictions = rf_model.predict(X_test)
        #train_predictions = rf_model.predict(X_test)
        #test_loss = np.mean((test_predictions - y_test)**2)  # Mean squared error
        #print(f'Test Loss: {test_loss}')

        f = open(gwas_result_name, 'w')
        data = []
        for col, score in zip(snp_ids, rf_model.feature_importances_):
            f.write(str(col) + ' ' + str(score) + '\n')
            data.append((col, score))
        # Convert the list of tuples into a DataFrame
        df = pd.DataFrame(data, columns=['SNP', 'PValue'])
        df[['Chr', 'ChrPos']] = df['SNP'].str.split(':', expand=True)
        return df


    def run_xgboost(self, snp_data, pheno_data, snp_ids, test_size, estimators, gwas_result_name, bed_gp, pheno_gp,predict ):
        print (predict)
        snp_data[np.isnan(snp_data)] = -1

        # Split data into training and testing sets
        X_train, X_test, y_train, y_test = train_test_split(snp_data, pheno_data, test_size=test_size, random_state=42)

        # # Standardize the input features
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train)
        X_test = scaler.transform(X_test)

        # Define and train a Random Forest model
        xgb_model = xgb.XGBRegressor(n_estimators=estimators, learning_rate=0.1, max_depth=3, random_state=42)
        #xgb_model = xgb.XGBRegressor(n_estimators=estimators, random_state=42)



        # Fit the model to the training data
        xgb_model.fit(X_train, y_train)

        if predict:
            print('Start GP')

            bed_data = bed_gp.read().iid
            pheno_data = pheno_gp.read().iid
            df_gp = pd.DataFrame(bed_data, columns=['ID1', 'BED_ID2'])
            predicted_values = xgb_model.predict(bed_gp.read().val)
            df_gp['Predicted_Value'] = predicted_values
            #print (df_gp)
            df_pheno = pd.DataFrame(pheno_data, columns=['ID1', 'Pheno_ID2'])
            df_pheno['Pheno_Value'] = pheno_gp.read().val
            merged_df = df_gp.merge(df_pheno, on='ID1', how='outer')
            merged_df.to_csv('gp_out.csv', sep=';', index=False)

            print (merged_df)
            #print (bed_gp.read().iid)

            #print (df_gp)


            #test_predictions = xgb_model.predict(X_test)
            #train_predictions = xgb_model.predict(X_train)
            #print (y_test, test_predictions)
            #print (y_train, train_predictions)
        # Evaluate the model on the test set
        #test_predictions = rf_model.predict(X_test)
        #test_loss = np.mean((test_predictions - y_test)**2)  # Mean squared error
        #print(f'Test Loss: {test_loss}')

        f = open(gwas_result_name, 'w')
        data = []
        for col, score in zip(snp_ids, xgb_model.feature_importances_):
            f.write(str(col) + ' ' + str(score) + '\n')
            data.append((col, score))
        # Convert the list of tuples into a DataFrame
        df = pd.DataFrame(data, columns=['SNP', 'PValue'])
        df[['Chr', 'ChrPos']] = df['SNP'].str.split(':', expand=True)
        return df

