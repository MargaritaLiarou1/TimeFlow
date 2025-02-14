import os
import torch
import numpy as np
import pandas as pd
from torch.utils.data import DataLoader

np.random.seed(1)
torch.manual_seed(1)

class Training:

    """
    Training class implements the training for density estimation models.
    """

    def __init__(self, session, max_patience, num_epochs, model, optimizer, training_loader, val_loader, D):

        self.session = session
        # Maximum patience for early stopping
        self.max_patience = max_patience
        # Number of epochs for training
        self.num_epochs = num_epochs
        self.model = model
        self.optimizer = optimizer
        # DataLoader for training data
        self.training_loader = training_loader
        # DataLoader for validation data
        self.val_loader = val_loader
        # Number of dimensions (CD markers/proteins)
        self.D = D
        # Use if test set is included, add it as argument
        # self.test_loader = test_loader

    def train(self):

        """
        Train the model using the specified parameters.

        Returns:
            np.ndarray: validation loss array over epochs
            int: epoch number needed for training to stop
        """

        # List to store validation loss values over epochs
        nll_val = []
        best_nll = 1000.0
        patience = 0

        # Main training loop
        for e in range(self.num_epochs):

            # Set the model to training mode
            self.model.train()
            # Iterate over batches in the training loader
            for indx_batch, batch in enumerate(self.training_loader):
                # Compute the loss for the current batch
                loss = self.model.forward(batch)

                # Back-propagate the loss and update the model parameters
                self.optimizer.zero_grad()
                loss.backward(retain_graph=True)
                self.optimizer.step()

            # Evaluate the model on the validation data
            loss_val = self.evaluation(self.val_loader, self.model, epoch=e)
            nll_val.append(loss_val)

            if loss_val < best_nll:
                print('Loss improved, new model is saved.')
                torch.save(self.model, f'{self.session}.model')
                best_nll = loss_val
                patience = 0
            else:
                patience += 1

            if patience > self.max_patience:
                break

        nll_val = np.asarray(nll_val)
        return nll_val, e

    def evaluation(self, val_loader, model, epoch):

        """
        Evaluate the model on the validation data.

        Args:
            val_loader: DataLoader for validation data
            model: model for evaluation
            epoch: current epoch number

        Returns:
            float: Validation loss
        """

        # Set the model to evaluation mode
        model.eval()
        loss_val = 0.0
        with torch.no_grad():
            # Iterate over batches in the validation loader
            for batch in val_loader:
                # Compute the loss for the current batch
                loss = model.forward(batch)
                loss_val += loss.item()

        # Return the average validation loss over all batches
        return loss_val / len(val_loader)

    def evaluate_models(models_directory, data_file_path, results_directory):

        """
        Evaluate all models in the specified directory and save the results.

        Args:
            models_directory: path to the directory containing model files.
            data_file_path: path to the data file.
            results_directory: directory to save evaluation results.
        """

        # Check the results directory exists
        if not os.path.exists(results_directory):
            os.makedirs(results_directory)

        # Iterate over model files in the directory
        for model_file in os.listdir(models_directory):
            if model_file.endswith(".model"):
                model_path = os.path.join(models_directory, model_file)
                model = torch.load(model_path)
                model.eval()

                # Load evaluation data
                flowcyteval = pd.read_csv(data_file_path)
                flowcyteval = flowcyteval.iloc[:, :self.D].to_numpy()
                flowcyteval = torch.tensor(flowcyteval[:, :].astype(np.float32))

                # Create DataLoader for evaluation data
                evalLoader = DataLoader(flowcyteval, batch_size=flowcyteval.shape[0], shuffle=False)

                # Evaluate model on the evaluation data
                for indx_batch, batch in enumerate(evalLoader):
                    output = model.forward(batch)

                # Convert model output to numpy array and create DataFrame
                output_array = output.detach().numpy()
                nf_pdf_df = pd.DataFrame(output_array, columns=["pdf"])

                # Save evaluation results to the results directory
                model_name = os.path.splitext(model_file)[0]
                save_path = os.path.join(results_directory, f"real_nvp_pdf_{model_name}.csv")
                nf_pdf_df.to_csv(save_path, index=False)
                print(f"Results for model: {model_name} are saved")

