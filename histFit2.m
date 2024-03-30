function histFit2(x,y,regions)

perfThresh = mean(regions(regions > 0));

lowPerfMask = regions < perfThresh;
highPerfMask = regions > perfThresh;
slopeLow = regions.*lowPerfMask;
slopeHigh = regions.*highPerfMask;

muLow = mean(slopeLow(slopeLow > 0));
sigmaLow = std(slopwLow(slopeLow > 0));
muHigh = mean(slopeHigh(slopeHigh >0));
sigmaHigh = std(slopeHigh(slopeHigh > 0));
X = x;
Y =y;

% Convert X and Y into a table, which is the form fitnlm() likes the input data to be in.
tbl = table(X', Y');
% Define the model as Y = a + b*x + c*exp(-(x-d)^2/e) + d * exp(-(x-f)^2/g)
% Note how this "x" of modelfun is related to big X and big Y.
% x(:, 1) is actually X and x(:, 2) is actually Y - the first and second columns of the table.
modelfun = @(b,x) b(1) + b(2) * x(:, 1) + b(3) * exp(-(x(:, 1) - b(4)).^2/b(5)) + b(6) * exp(-(x(:, 1) - b(7)).^2/b(8));  
beta0 = [6, 0.1, 35, 10, 3, 25, 50, 9]; % Guess values to start with.  Just make your best guess.
% Now the next line is where the actual model computation is done.
mdl = fitnlm(tbl, modelfun, beta0);





end