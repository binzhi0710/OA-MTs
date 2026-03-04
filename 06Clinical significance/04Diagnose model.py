import pandas as pd
import numpy as np
from sklearn.model_selection import KFold
from sklearn.metrics import roc_auc_score, roc_curve, recall_score, precision_score, f1_score
import matplotlib.pyplot as plt
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression, RidgeClassifierCV
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.neural_network import MLPClassifier

# Set global random seed
SEED = 2024
np.random.seed(SEED)

# Read feature genes
with open("feature_genes.txt") as f:
    feature_genes = [x.strip() for x in f.readlines()]

# Data loading function
def load_data(expr_file, group_file):
    expr = pd.read_csv(expr_file, sep='\t', index_col=0).T
    valid_genes = [g for g in feature_genes if g in expr.columns]
    expr = expr[valid_genes].fillna(expr[valid_genes].mean())

    groups = pd.read_excel(group_file)[['name', 'type']].set_index('name')
    groups['type'] = groups['type'].map({'control': 0, 'OA': 1}).astype(int)

    data = expr.join(groups, how='inner')
    return data[valid_genes], data['type']

# Load data
X_train, y_train = load_data("trainingsample_z_score_normalized_data.txt", "group_info.xlsx")
X_valA, y_valA = load_data("GSE55235.txt", "GSE55235_group_info.xlsx")
X_valB, y_valB = load_data("GSE114007.txt", "GSE114007_group_info.xlsx")

# Initialize models (with random seed)
models = {
    "Naive Bayes": GaussianNB(),
    "Logistic Regression": LogisticRegression(max_iter=1000, random_state=SEED),
    "Lasso": LogisticRegression(penalty='l1', solver='liblinear', max_iter=1000, random_state=SEED),
    "RidgeCV": RidgeClassifierCV(),
    "ElasticNet": LogisticRegression(penalty='elasticnet', solver='saga', l1_ratio=0.5, max_iter=1000, random_state=SEED),
    "LDA": LinearDiscriminantAnalysis(),
    "MLP": MLPClassifier(max_iter=1000, hidden_layer_sizes=(100,), random_state=SEED)
}

# Store results
results = []
train_predictions = {}

# Train and evaluate models
for name, model in models.items():
    print(f"Training {name}...")

    # Five-fold cross-validation
    kf = KFold(n_splits=5, shuffle=True, random_state=SEED)
    all_y_test = []
    all_y_proba = []
    all_y_pred = []

    for train_idx, test_idx in kf.split(X_train):
        X_tr, X_te = X_train.iloc[train_idx], X_train.iloc[test_idx]
        y_tr, y_te = y_train.iloc[train_idx], y_train.iloc[test_idx]

        model.fit(X_tr, y_tr)

        # Get predictions
        if hasattr(model, 'predict_proba'):
            y_proba = model.predict_proba(X_te)[:, 1]
        else:
            y_proba = model.decision_function(X_te)
        y_pred = model.predict(X_te)

        all_y_test.extend(y_te)
        all_y_proba.extend(y_proba)
        all_y_pred.extend(y_pred)

    # Calculate cross-validation metrics
    cv_metrics = {
        'AUC': roc_auc_score(all_y_test, all_y_proba),
        'Recall': recall_score(all_y_test, all_y_pred),
        'Precision': precision_score(all_y_test, all_y_pred),
        'F1': f1_score(all_y_test, all_y_pred)
    }

    # Store training set predictions for plotting
    train_predictions[name] = (all_y_test, all_y_proba)

    # Full training
    model.fit(X_train, y_train)

    # Evaluation function
    def evaluate_model(X, y):
        if hasattr(model, 'predict_proba'):
            y_proba = model.predict_proba(X)[:, 1]
        else:
            y_proba = model.decision_function(X)
        y_pred = model.predict(X)
        return {
            'AUC': roc_auc_score(y, y_proba),
            'Recall': recall_score(y, y_pred),
            'Precision': precision_score(y, y_pred),
            'F1': f1_score(y, y_pred)
        }

    # Record results for validation sets
    valA_metrics = evaluate_model(X_valA, y_valA)
    valB_metrics = evaluate_model(X_valB, y_valB)
    results.append({
        'Model': name,
        'Train AUC': cv_metrics['AUC'],
        'Train Recall': cv_metrics['Recall'],
        'Train Precision': cv_metrics['Precision'],
        'Train F1': cv_metrics['F1'],
        'ValA AUC': valA_metrics['AUC'],
        'ValA Recall': valA_metrics['Recall'],
        'ValA Precision': valA_metrics['Precision'],
        'ValA F1': valA_metrics['F1'],
        'ValB AUC': valB_metrics['AUC'],
        'ValB Recall': valB_metrics['Recall'],
        'ValB Precision': valB_metrics['Precision'],
        'ValB F1': valB_metrics['F1']
    })

# Save results
pd.DataFrame(results).to_excel("model_performance_v2.xlsx", index=False)

# Plot ROC curves
def plot_roc(dataset_name, models_predictions, models, X_data=None, y_data=None):
    plt.figure(figsize=(8, 8))
    for name in models:
        if dataset_name == "Training Set":
            y_true, y_score = models_predictions[name]
        else:  # Validation Set
            y_true = y_data
            y_score = models[name].predict_proba(X_data)[:, 1] if hasattr(models[name], 'predict_proba') else models[name].decision_function(X_data)

        fpr, tpr, _ = roc_curve(y_true, y_score)
        auc = roc_auc_score(y_true, y_score)
        plt.plot(fpr, tpr, label=f'{name} (AUC={auc:.2f})')

    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlabel('1-Specificity')
    plt.ylabel('Sensitivity')
    plt.title(f'ROC Curves ({dataset_name})')
    plt.legend(loc='lower right')
    plt.savefig(f'{dataset_name.replace(" ", "_")}_ROC.pdf', bbox_inches='tight')
    plt.close()

# Plot ROC curves for training set and both validation sets
plot_roc("Training Set", train_predictions, models)
plot_roc("Validation Set A", None, models, X_valA, y_valA)
plot_roc("Validation Set B", None, models, X_valB, y_valB)

print("All models trained! Results saved to model_performance_v2.xlsx")