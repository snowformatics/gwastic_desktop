import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import xgboost as xgb
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.metrics import accuracy_score, mean_squared_error


class GWASAI:
    def run_random_forest(self, snp_data, pheno_data, snp_ids, test_size, estimators):

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
        #test_loss = np.mean((test_predictions - y_test)**2)  # Mean squared error
        #print(f'Test Loss: {test_loss}')

        f = open('out_rf_test.txt', 'w')
        data = []
        for col, score in zip(snp_ids, rf_model.feature_importances_):
            f.write(str(col) + ' ' + str(score) + '\n')
            data.append((col, score))
        # Convert the list of tuples into a DataFrame
        df = pd.DataFrame(data, columns=['SNP', 'Value'])
        df[['Chr', 'ChrPos']] = df['SNP'].str.split(':', expand=True)
        return df


def plot():
    #df = pd.read_csv('out_xgboos_bridge2.txt', delimiter=' ', header=None)
    df = pd.read_csv('C:/Users/lueck/PycharmProjects/genbank2/test/_out_nn_bridge100.txt', delimiter=' ', header=None)

    df.columns = ['snp', 'value']
    #print (df.dtypes)
    #print(df)
    #df = df.sort_values(by=['value'], ascending=False)
    #ax1 = df.plot.scatter(x='value',
                        #  y='snp')
    #plt.savefig('out_xgboos_bridge.png')
    feature_list = df['value'].to_numpy()
    plt.plot(feature_list)
    plt.show()
    #print (feature_list)
    #plt.axhline(y=feature_list, color='b', linestyle='--')
    # plt.ylabel('saliency value', fontdict=None, labelpad=None, fontsize=15)
    #
    # plt.xlabel('SNPs', fontdict=None, labelpad=None, fontsize=15)
    #
    #plt.savefig('sal.png')
    #

    # print (len(xgb_model.feature_importances_))
    # xgb.plot_importance(xgb_model, max_num_features=20)
    # plt.rcParams['figure.figsize'] = [50, 50]
    # plt.savefig('rf.png')

#run_tree()
#plot()
