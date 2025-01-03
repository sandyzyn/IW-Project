import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.decomposition import PCA


file_path = 'orco gene - All comparison.csv'  
data = pd.read_csv(file_path)

#break down data
all_exons = []
num_exons = 8

for i in range(num_exons):
    start_col = i * 6  
    exon_data = data.iloc[:, [start_col + 1, start_col + 2, start_col + 3, start_col + 4]].dropna()
    exon_data.columns = ['Benchling_MIT', 'Geneious', 'Guidescan2', 'Benchling_CFD']

    for col in ['Benchling_MIT', 'Geneious', 'Guidescan2', 'Benchling_CFD']:
        exon_data[col] = pd.to_numeric(exon_data[col], errors='coerce')
    exon_data = exon_data.dropna()

    exon_data['Exon'] = f'Exon {i + 1}'
    all_exons.append(exon_data)

combined_data = pd.concat(all_exons, ignore_index=True)

X = combined_data[['Benchling_MIT', 'Geneious', 'Guidescan2', 'Benchling_CFD']]

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

reduced = PCA(n_components=1).fit_transform(X_scaled)

# 0-100 range
min_max_scaler = MinMaxScaler(feature_range=(0, 100))
y = min_max_scaler.fit_transform(reduced)

def train(X, y):
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    model = LinearRegression()
    model.fit(X_train, y_train.ravel())
    y_hat = model.predict(X_test)
    mse = mean_squared_error(y_test, y_hat)
    r2 = r2_score(y_test, y_hat)
    return model, mse, r2, X_test, y_test, y_hat

def predict(model):
    try:
        print("\nEnter new scores to predict the unified score:")
        benchling_mit = float(input("Benchling_MIT: "))
        geneious = float(input("Geneious: "))
        guidescan2 = float(input("Guidescan2: "))
        benchling_cfd = float(input("Benchling_CFD: "))
        input_data = pd.DataFrame([[benchling_mit, geneious, guidescan2, benchling_cfd]],
                                  columns=['Benchling_MIT', 'Geneious', 'Guidescan2', 'Benchling_CFD'])
        predicted_score = model.predict(input_data)[0]
        print(f"\nPredicted Unified Score: {predicted_score}")
    except ValueError:
        print("Invalid input. Please enter numeric values.")

def feature_importance(model, X):
    importance = model.coef_  
    features = X.columns
    plt.barh(features, importance)
    plt.xlabel("Importance")
    plt.ylabel("Software")
    plt.title("Feature Importance in Predicting Unified Score")
    plt.show()

if __name__ == "__main__":
    model, mse, r2, X_test, y_test, y_hat = train(X, y)
    print(f"Mean Squared Error: {mse}")
    print(f"R² Score: {r2}")

    while True:
        print("\nOptions:")
        print("1. Predict Score")
        print("2. Show Feature Importance")
        print("3. Exit")

        choice = input("Choose an option: ").strip()
        if choice == "1":
            predict(model)
        elif choice == "2":
            feature_importance(model, X)
        elif choice == "3":
            print("Exiting program.")
            break
        else:
            print("Invalid choice. Please try again.")
