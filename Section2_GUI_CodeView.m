classdef CIE327_PROJECT < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        VarianceofYLabel                matlab.ui.control.Label
        VarianceofXLabel                matlab.ui.control.Label
        MeanofYLabel                    matlab.ui.control.Label
        MeanofXLabel                    matlab.ui.control.Label
        FileNameWithExtensionEditField  matlab.ui.control.EditField
        CalculateCovarianceandCorrelationCoefficientButton  matlab.ui.control.Button
        CalculateandPlotMarginalDistributionButton  matlab.ui.control.Button
        CalculateandPlotJointDistributionButton  matlab.ui.control.Button
        FileNameWithExtensionEditFieldLabel  matlab.ui.control.Label
        CorrelationCoefficientLabel     matlab.ui.control.Label
        CovariancexyLabel               matlab.ui.control.Label
        UIAxes3                         matlab.ui.control.UIAxes
        UIAxes2                         matlab.ui.control.UIAxes
        UIAxes                          matlab.ui.control.UIAxes
    end

    
    properties (Access = private)
        Property % Description
        X % Data for X
        Y % Data for Y
        xEdges % Region edges for X
        yEdges % Region edges for Y
        jointProbabilities % Stores the joint probability distribution matrix
        xCenters % Region centers for X
        yCenters % Region centers for Y
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: CalculateandPlotJointDistributionButton
        function CalculateandPlotJointDistributionButtonPushed(app, event)

            % Ensure X and Y are available 
            X = app.X;
            Y = app.Y;

            % Ensure X and Y are not empty
            if isempty(X) || isempty(Y)
                uialert(app.UIFigure, 'Data is not loaded. Please load a valid file.', 'Error');
                return;
            end

            % Ensure the region edges have at least 2 values
            if numel(app.xEdges) < 2 || numel(app.yEdges) < 2
                error('xEdges or yEdges must have at least 2 values.');
            end

            % Sort region edges (already defined when loading the .mat file)
            xEdges = sort(app.xEdges);
            yEdges = sort(app.yEdges);

            % Check the size of the region edges and ranges (for debugging)
            disp(['Number of xEdges: ', num2str(numel(xEdges))]);
            disp(['Number of yEdges: ', num2str(numel(yEdges))]);
            disp(['X range: ', num2str(min(X)), ' to ', num2str(max(X))]);
            disp(['Y range: ', num2str(min(Y)), ' to ', num2str(max(Y))]);

            % Descritize X and Y into indices 
            xRegionIndices = discretize(X, xEdges);
            yRegionIndices = discretize(Y, yEdges);

            % Calculate the Joint Frequency Distribution
            [jointCounts, ~, ~] = histcounts2(X, Y, xEdges, yEdges);

            % Check jointCounts (for debugging)
            disp('Joint Counts:');
            disp(jointCounts);

            if sum(jointCounts(:)) == 0
                uialert(app.UIFigure, 'No data points fall into the specified bins.', 'Error');
                return;
            end
            
            % Normalize to obtain the joint probability distribution 
            jointProbabilities = jointCounts / sum(jointCounts(:));

            % Check if jointProbabilities sum to 1
             disp(['Sum of joint probabilities: ', num2str(sum(jointProbabilities(:)))]);

             % Generate the meshgrid for plotting
             xCenters = (xEdges(1:end-1) + xEdges(2:end)) / 2;
             yCenters = (yEdges(1:end-1) + yEdges(2:end)) / 2;
            [XGrid, YGrid] = meshgrid(xCenters, yCenters);

            % Ensure jointProbabilities has the correct dimensions
            [nRows, nCols] = size(jointCounts);
            if numel(xEdges) - 1 ~= nCols || numel(yEdges) - 1 ~= nRows
                error('The dimensions of the joint probabilities do not match the number of the regions.');
            end    

            % Save joint probabilities and region centers to app properties
            app.jointProbabilities = jointProbabilities; 
            app.xCenters = (xEdges(1:end-1) + xEdges(2:end)) / 2; 
            app.yCenters = (yEdges(1:end-1) + yEdges(2:end)) / 2; 


            % Plot the joint probability distribution 
            axes(app.UIAxes);
            surf(app.UIAxes, XGrid, YGrid, jointProbabilities, 'EdgeColor', 'none');
            colorbar(app.UIAxes);
            xlabel(app.UIAxes, 'X regions');
            ylabel(app.UIAxes, 'Y regions');
            zlabel(app.UIAxes, 'Joint Probability');
            title(app.UIAxes, 'Joint Probability Distribution');

            % Adjust the view angle for better visualization
            view(app.UIAxes, [120 30]);

        end

        % Button pushed function: 
        % CalculateandPlotMarginalDistributionButton
        function CalculateandPlotMarginalDistributionButtonPushed(app, event)
            % Ensure jointProbabilities is available
            if isempty(app.jointProbabilities)
                uialert(app.UIFigure, 'Joint probabilities are not calculated. Please calculate them first.', 'Error');
                return;
            end    

            % Ensure X and Y centers are available
            if isempty(app.xCenters) || isempty(app.yCenters)
                uialert(app.UIFigure, 'X or Y centers are not available. Please check your data.', 'Error');
                return;
            end;    

            % Calculate marginal distributions
            marginalX = sum(app.jointProbabilities, 2); 
            marginalY = sum(app.jointProbabilities, 1); 

            % Normalize marginals to ensure they sum to 1
            marginalX = marginalX / sum(marginalX);
            marginalY = marginalY / sum(marginalY); 

            % Plot Marginal Distribution of X
            axes(app.UIAxes2); 
            plot(app.UIAxes2, app.xCenters, marginalX, 'LineWidth', 2, 'Color', [0.2, 0.6, 0.8]);
            xlabel(app.UIAxes2, 'X'); 
            ylabel(app.UIAxes2, 'Probability'); 
            title(app.UIAxes2, 'Marginal Distribution of X'); 
            grid(app.UIAxes2, 'on'); 

            % Plot Marginal Distribution of Y
            axes(app.UIAxes3); 
            plot(app.UIAxes3, app.yCenters, marginalY, 'LineWidth', 2, 'Color', [0.8, 0.4, 0.2]);
            xlabel(app.UIAxes3, 'Y'); 
            ylabel(app.UIAxes3, 'Probability'); 
            title(app.UIAxes3, 'Marginal Distribution of Y'); 
            grid(app.UIAxes3, 'on');

        end

        % Button pushed function: 
        % CalculateCovarianceandCorrelationCoefficientButton
        function CalculateCovarianceandCorrelationCoefficientButtonPushed(app, event)

            % Ensure jointProbabilities are available
            if isempty(app.jointProbabilities)
                uialert(app.UIFigure, 'Joint probabilities are not calculated. Please calculate them first.', 'Error');
                return;
            end

            % Ensure X and Y centers are available
            if isempty(app.xCenters) || isempty(app.yCenters)
                uialert(app.UIFigure, 'Region centers are not available. Please check your data.', 'Error');
                return;
            end

            % Normalize joint probabilities to ensure they sum to 
            totalProbability = sum(app.jointProbabilities, 'all');
            if abs(totalProbability - 1) > 1e-6
                app.jointProbabilities = app.jointProbabilities / totalProbability;
            end

            % Marginal distributions
            marginalX = sum(app.jointProbabilities, 2); 
            marginalY = sum(app.jointProbabilities, 1)';

            % Ensure vectors are aligned
            xCenters = app.xCenters(:); 
            yCenters = app.yCenters(:); 

            % Compute expected values (means)
            meanX = sum(xCenters .* marginalX); % E[X] = sum(x * f(x))
            meanY = sum(yCenters .* marginalY); % E[Y] = sum(y * f(y))

            % Compute expected values of x^2 and y^2
            expectedX2 = sum((xCenters.^2) .* marginalX); % E[X^2] = sum(x^2 * f(x))
            expectedY2 = sum((yCenters.^2) .* marginalY); % E[Y^2] = sum(y^2 * f(y))

            % Compute variances
            varianceX = expectedX2 - meanX^2; % Variance of X = E[X^2] - (E[X])^2
            varianceY = expectedY2 - meanY^2; % Variance of Y = E[Y^2] - (E[Y])^2

           % Validate variances to avoid division by zero
           if varianceX <= 0 || varianceY <= 0
               uialert(app.UIFigure, 'Variance cannot be zero or negative. Check your data.', 'Error');
               return;
           end

           % Compute expected the expected value of xy
           expectedXY = 0; 
           for i = 1:numel(xCenters)
               for j = 1:numel(yCenters)
                   expectedXY = expectedXY + (xCenters(i) * yCenters(j) * app.jointProbabilities(i, j));
               end
           end

           % Compute covariance
           covarianceXY = expectedXY - (meanX * meanY); % Covariance = E[XY] - E[X] * E[Y]

           % Compute correlation coefficient
           correlationCoefficient = covarianceXY / sqrt(varianceX * varianceY);

           % Validate correlation coefficient
           if abs(correlationCoefficient) > 1
               uialert(app.UIFigure, 'Correlation coefficient exceeds valid range. Check your data.', 'Error');
               return;
           end

           % Update UI labels with results
           app.MeanofXLabel.Text = sprintf('Mean of X: %.4f', meanX);
           app.MeanofYLabel.Text = sprintf('Mean of Y: %.4f', meanY);
           app.VarianceofXLabel.Text = sprintf('Variance of X: %.4f', varianceX);
           app.VarianceofYLabel.Text = sprintf('Variance of Y: %.4f', varianceY);
           app.CovariancexyLabel.Text = sprintf('Covariance(xy): %.4f', covarianceXY);
           app.CorrelationCoefficientLabel.Text = sprintf('Correlation Coefficient: %.4f', correlationCoefficient);

        end

        % Value changed function: FileNameWithExtensionEditField
        function FileNameWithExtensionEditFieldValueChanged(app, event)
            % Get the file name entered by the user
            filename = app.FileNameWithExtensionEditField.Value;
      try
            % Load the file
            data = load(filename);

            % Ensure the loaded data contains the variable 'XY'
            if isfield(data, 'XY')
                % Extract the XY matrix from the loaded data structure 
                XY = data.XY;

                % Extract X and Y from the matrix
                X = XY(1, :);
                Y = XY(2, :);

                % Store X and Y in app properties for later use
                app.X = X;
                app.Y = Y;

                % Define the ranges of X and Y
                xMin = min(X);
                xMax = max(X);
                yMin = min(Y);
                yMax = max(Y);

                % Define the number of regions 
                numRegions = 500;

                % Create region edges and store for later use
                app.xEdges = linspace(xMin, xMax, numRegions + 1);
                app.yEdges = linspace(yMin, yMax, numRegions + 1);

                % Notify the user of sucessful file entry 
                uialert(app.UIFigure, 'File loaded successfully!','Success');
            else 
                % Error handling for missing variable 'XY'
                uialert(app.UIFigure,'Error, the file does not contain the variable XY.','Error');
            end
      catch
            % Error handling
            uialert(app.UIFigure, 'Error, Please check the filename and format.','Error');
      end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [6 6 1000 750];
            app.UIFigure.Name = 'MATLAB App';

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            title(app.UIAxes, 'Joint Probability Distribution')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.Position = [20 340 460 220];

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.UIFigure);
            title(app.UIAxes2, 'Marginal Distribution of X')
            xlabel(app.UIAxes2, 'X')
            ylabel(app.UIAxes2, 'Y')
            zlabel(app.UIAxes2, 'Z')
            app.UIAxes2.Position = [20 110 460 220];

            % Create UIAxes3
            app.UIAxes3 = uiaxes(app.UIFigure);
            title(app.UIAxes3, 'Marginal Distribution of Y')
            xlabel(app.UIAxes3, 'X')
            ylabel(app.UIAxes3, 'Y')
            zlabel(app.UIAxes3, 'Z')
            app.UIAxes3.Position = [500 110 460 220];

            % Create CovariancexyLabel
            app.CovariancexyLabel = uilabel(app.UIFigure);
            app.CovariancexyLabel.FontSize = 13;
            app.CovariancexyLabel.Position = [20 70 200 22];
            app.CovariancexyLabel.Text = 'Covariance(xy): ';

            % Create CorrelationCoefficientLabel
            app.CorrelationCoefficientLabel = uilabel(app.UIFigure);
            app.CorrelationCoefficientLabel.FontSize = 13;
            app.CorrelationCoefficientLabel.Position = [240 70 200 22];
            app.CorrelationCoefficientLabel.Text = 'Correlation Coefficient:';

            % Create FileNameWithExtensionEditFieldLabel
            app.FileNameWithExtensionEditFieldLabel = uilabel(app.UIFigure);
            app.FileNameWithExtensionEditFieldLabel.HorizontalAlignment = 'right';
            app.FileNameWithExtensionEditFieldLabel.FontSize = 13;
            app.FileNameWithExtensionEditFieldLabel.Position = [20 669 161 22];
            app.FileNameWithExtensionEditFieldLabel.Text = 'File Name (With Extension)';

            % Create CalculateandPlotJointDistributionButton
            app.CalculateandPlotJointDistributionButton = uibutton(app.UIFigure, 'push');
            app.CalculateandPlotJointDistributionButton.ButtonPushedFcn = createCallbackFcn(app, @CalculateandPlotJointDistributionButtonPushed, true);
            app.CalculateandPlotJointDistributionButton.FontSize = 13;
            app.CalculateandPlotJointDistributionButton.Position = [20 590 280 40];
            app.CalculateandPlotJointDistributionButton.Text = 'Calculate and Plot Joint Distribution';

            % Create CalculateandPlotMarginalDistributionButton
            app.CalculateandPlotMarginalDistributionButton = uibutton(app.UIFigure, 'push');
            app.CalculateandPlotMarginalDistributionButton.ButtonPushedFcn = createCallbackFcn(app, @CalculateandPlotMarginalDistributionButtonPushed, true);
            app.CalculateandPlotMarginalDistributionButton.FontSize = 13;
            app.CalculateandPlotMarginalDistributionButton.Position = [361 590 280 40];
            app.CalculateandPlotMarginalDistributionButton.Text = 'Calculate and Plot Marginal Distribution';

            % Create CalculateCovarianceandCorrelationCoefficientButton
            app.CalculateCovarianceandCorrelationCoefficientButton = uibutton(app.UIFigure, 'push');
            app.CalculateCovarianceandCorrelationCoefficientButton.ButtonPushedFcn = createCallbackFcn(app, @CalculateCovarianceandCorrelationCoefficientButtonPushed, true);
            app.CalculateCovarianceandCorrelationCoefficientButton.FontSize = 13;
            app.CalculateCovarianceandCorrelationCoefficientButton.Position = [680 590 302 40];
            app.CalculateCovarianceandCorrelationCoefficientButton.Text = 'Calculate Covariance and Correlation Coefficient';

            % Create FileNameWithExtensionEditField
            app.FileNameWithExtensionEditField = uieditfield(app.UIFigure, 'text');
            app.FileNameWithExtensionEditField.ValueChangedFcn = createCallbackFcn(app, @FileNameWithExtensionEditFieldValueChanged, true);
            app.FileNameWithExtensionEditField.FontSize = 13;
            app.FileNameWithExtensionEditField.Position = [200 665 250 30];

            % Create MeanofXLabel
            app.MeanofXLabel = uilabel(app.UIFigure);
            app.MeanofXLabel.FontSize = 13;
            app.MeanofXLabel.Position = [548 517 180 22];
            app.MeanofXLabel.Text = 'Mean of X:';

            % Create MeanofYLabel
            app.MeanofYLabel = uilabel(app.UIFigure);
            app.MeanofYLabel.FontSize = 13;
            app.MeanofYLabel.Position = [761 517 221 22];
            app.MeanofYLabel.Text = 'Mean of Y: ';

            % Create VarianceofXLabel
            app.VarianceofXLabel = uilabel(app.UIFigure);
            app.VarianceofXLabel.FontSize = 13;
            app.VarianceofXLabel.Position = [548 393 180 22];
            app.VarianceofXLabel.Text = 'Variance of X:';

            % Create VarianceofYLabel
            app.VarianceofYLabel = uilabel(app.UIFigure);
            app.VarianceofYLabel.Position = [761 393 199 22];
            app.VarianceofYLabel.Text = 'Variance of Y:';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = CIE327_PROJECT

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end