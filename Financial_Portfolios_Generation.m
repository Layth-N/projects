clear
clc

% In this part we are calculating the sample mean and sample variance of
% our historical data.
% For the different portfolios calculations skip to the next part.

%define the series of CSV file names
Assets = {'F.csv', 'CAT.csv', 'DIS.csv', 'MCD.csv', 'KO.csv', 'PEP.csv', 'WMT.csv', 'C.csv', 'WFC.csv', 'JPM.csv', 'AAPL.csv', 'IBM.csv', 'PFE.csv', 'JNJ.csv', 'XOM.csv', 'MRO.csv', 'ED.csv', 'T.csv', 'VZ.csv', 'NEM.csv'};
Assets_count = numel(Assets);

%define a cell array to store the data from each CSV file
Assets_value = cell(1, Assets_count);
%loop through each CSV file
for i = 1:Assets_count
    %read the CSV file
    Assets_value{i} = readtable(Assets{i});
    %convert the first column to datetime format
    Assets_value{i}.Date = datetime(Assets_value{i}.Date);
    %extract month and year information
    Assets_value{i}.Month = month(Assets_value{i}.Date);
end

%define a cell array to store return of each end of month
end_of_month = cell(1, Assets_count);
%getting total number of months
monthYear = year(Assets_value{1}.Date) * 100 + month(Assets_value{1}.Date);
total_months_count = size(histcounts(monthYear, unique(monthYear)),2);
zeros_1 = zeros(total_months_count,1);

%getting the end of month values
for i = 1:Assets_count
    past_month = 12;
    p = 1;
    end_of_month{i} = table(zeros_1, zeros_1, 'VariableNames', {'End_of_Month', 'Return'});
    for j = 1:height(Assets_value{i})
        current_month = Assets_value{i}.Month(j);
        if current_month ~= past_month
            end_of_month{i}.End_of_Month(p) = past_month;
            end_of_month{i}.Return(p) = Assets_value{i}.AdjClose(j-1);
            p = p+1;
        end
        past_month = Assets_value{i}.Month(j);
    end
end

%define a cell array to store months returns
months_returns = cell(1, Assets_count);
%define a matrix to store arithmatic average
Average_returns = zeros(1, Assets_count);
%define a matrix to store Assets returns
Assets_returns = zeros(1, Assets_count);

%Getting the rolling span value
Rolling_Span = 21;
June_2011_Position = 43;
zeros_2 = zeros(Rolling_Span,1);

for i = 1:Assets_count
    months_returns{i} = table(zeros_2, zeros_2, zeros_2, 'VariableNames', {'Month', 'Return', 'Return_minus_Average'});
    for j = June_2011_Position-Rolling_Span+1:June_2011_Position
        months_returns{i}.Month(j-June_2011_Position+Rolling_Span) = end_of_month{i}.End_of_Month(j);
        %calculating monthly_returns(i) = (last_day_of_month(i)/ last_day_of_month(i-1))-1
        months_returns{i}.Return(j-June_2011_Position+Rolling_Span) = ((end_of_month{i}.Return(j)/end_of_month{i}.Return(j-1))-1);
    end
    %calcualting arithmatic average = avg(r_it)
    Average_returns(i) = mean(months_returns{i}.Return);
    %calcualting geometric expected return of asset i = ((product(1+r_it))^(1/T))-1
    Assets_returns(i) = ((prod((months_returns{i}.Return)+1))^(1/Rolling_Span))-1;
    
    %calculating r_it - avg(r_it) (for covariance calculations later)
    for k = 1:Rolling_Span
        months_returns{i}.Return_minus_Average(k) = months_returns{i}.Return(k)-Average_returns(i);
    end
end

%define a matrix to store Assets covariance
Assets_covariance = zeros(Assets_count, Assets_count);
for i = 1:Assets_count
    for j = 1:Assets_count
        Difference_product = 1;
        Difference_product_sum = 0;
        if Assets_covariance(i,j) == 0 %to avoid calculating the covariance twice
            for k = 1:Rolling_Span
                %[r_it - avg(r_it)] * [r_jt - avg(r_jt)]
                Difference_product = months_returns{i}.Return_minus_Average(k) * months_returns{j}.Return_minus_Average(k);
                %sum{[r_it - avg(r_it)] * [r_jt - avg(r_jt)]}
                Difference_product_sum = Difference_product_sum + Difference_product;
            end
            %covariance between assets i and j = sum{[r_it - avg(r_it)] * [r_jt - avg(r_jt)]}/T
            Assets_covariance(i,j) = Difference_product_sum/Rolling_Span;
            Assets_covariance(j,i) = Assets_covariance(i,j);
        end
    end
end
  
% Compute the eigenvalues of the matrix
eigenvalues = eig(Assets_covariance);
   
% Check if all eigenvalues are non-negative
if all(eigenvalues >= 0)
    disp('The matrix is PSD.');
else
    disp('The matrix is not PSD.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculating Lambda
MarketCap_table = readtable('Market Cap.csv');

Market_Weights = zeros(1, Assets_count);
Market_returns = zeros(1, Assets_count);

for i = 1:Assets_count
    Market_Weights(i) = MarketCap_table.Average(i)/sum(MarketCap_table.Average); %assets weights = asset market cap/total market cap
    Market_returns(i) = Market_Weights(i) * Average_returns(i);
end

r_mkt = sum(Market_returns); %market return = weighted sum of each asset return

Market_Variance = 0;

for i = 1:Assets_count
    for j = 1:Assets_count
        Market_Variance = Market_Variance + Market_Weights(i) * Market_Weights(j) * Assets_covariance(i,j); %Market Variance = Sum(wi wj segam ij)
    end
end

r_f = 0.04/12;

lambda = (r_mkt - r_f) / Market_Variance;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%MVO

n=Assets_count;
mu=Average_returns; % expected returns of assets

Q= lambda * Assets_covariance; % covariance matrix
c= -0.5 * mu; %linear coefficients
Aeq =[ones(1,n)]; %equal A matrix
beq =[1]; %equal b vector
ub = inf * ones(n,1); %upper bound
lb = -inf * ones(n,1); %lower bound with short selling

%Risk Adjusted MVO
[MVO_X(1,:),fval_MVO(1,1)] = quadprog (Q, c, [], [], Aeq, beq, lb, ub);
MVO_r = MVO_X(1,:)*mu';
MVO_VAR = MVO_X(1,:) * Q * MVO_X(1,:)';
MVO_SR = (MVO_r - r_f)/((MVO_VAR)^0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Robust MVO

%Calculating delta
Assets_std = zeros(1, Assets_count);
Theta_sqr = zeros(1, Assets_count);
period_sqr = Rolling_Span^0.5;

for i = 1:Assets_count
    Assets_std(i) = Assets_covariance(i,i)^0.5;
    Theta_sqr(i) = Assets_std(i)/period_sqr;
end

delta_90 = 1.645*Theta_sqr;
delta_95 = 1.96*Theta_sqr;

%Adjusting the Covariance Matrix
Q_Robust = zeros(2*Assets_count, 2*Assets_count);
for i = 1:Assets_count
    for j = 1:Assets_count
        Q_Robust(i,j) = Assets_covariance(i,j);
    end
end

for i = Assets_count+1:2*Assets_count
    for j = Assets_count+1:2*Assets_count
        Q_Robust(i,j) = 0;
    end
end

%Ajusting the Unequal Constraint Matrix
A = zeros(2*Assets_count, 2*Assets_count);

for i = 1:Assets_count
    A(i, i) = 1;
    A(i, Assets_count+i) = -1;
    A(Assets_count+i, i) = -1;
    A(Assets_count+i, Assets_count+i) = -1;
end

Q_Robust = lambda * Q_Robust; % covariance matrix
c_90= [-0.5 * mu , 0.5 * delta_90]; %linear coefficients
c_95= [-0.5 * mu , 0.5 * delta_95]; %linear coefficients
A = A;
b = [zeros(1,2*n)];
Aeq_Robust =[ones(1,n) , zeros(1,n)]; %equal A matrix
beq =[1]; %equal b vector
ub_Robust = inf * ones(2*n,1); %upper bound
lb_Robust = -inf * ones(2*n,1); %lower bound with short selling

%Robust Risk Adjusted MVO 90%
[MVO_R_90_X(1,:),fval_MVO_R_90(1,1)] = quadprog (Q_Robust, c_90, A, b, Aeq_Robust, beq, lb_Robust, ub_Robust);
MVO_R_90_r = MVO_R_90_X(1:n)* mu';
MVO_R_90_VAR = MVO_R_90_X(1:n) * Q * MVO_R_90_X(1:n)';
MVO_R_90_SR = (MVO_R_90_r - r_f)/((MVO_R_90_VAR)^0.5);

%Robust Risk Adjusted MVO 95%
[MVO_R_95_X(1,:),fval_MVO_R_95(1,1)] = quadprog (Q_Robust, c_95, A, b, Aeq_Robust, beq, lb_Robust, ub_Robust);
MVO_R_95_r = MVO_R_95_X(1:n)* mu';
MVO_R_95_VAR = MVO_R_95_X(1:n) * Q * MVO_R_95_X(1:n)';
MVO_R_95_SR = (MVO_R_95_r - r_f)/((MVO_R_95_VAR)^0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Risk Parity

%Storing the Sparse Semmetric Matrix
R = cell(1, Assets_count);

for i = 1:Assets_count
    for j = 1:Assets_count
        R{i}(i,j) = 0.5*Assets_covariance(i,j);
        R{i}(j,i) = 0.5*Assets_covariance(i,j);
        if i == j
            R{i}(i,j) = Assets_covariance(i,j);
        end
    end
end

x0 = [zeros(n, 1); 0.1];
Aeq_RP =[ones(1,n),0]; %equal A matrix
lb_RP = zeros(n, 1); %no short selling lower bound
objfun = @(vars) objective(vars, R, n);
options = optimoptions('fmincon', 'Display', 'iter');
[Risk_Parity_X, fval_Risk_Parity, exitflag, output] = fmincon(objfun, x0, [], [], Aeq_RP, beq, lb_RP, [], [], options);

Risk_Parity_r = Risk_Parity_X(1:n)'* mu';
Risk_Parity_VAR = Risk_Parity_X(1:n)' * Q * Risk_Parity_X(1:n);
Risk_Parity_SR = (Risk_Parity_r - r_f)/((Risk_Parity_VAR)^0.5);

%Checking the risk parity solution
RP_CHECK = zeros(1, n);
for i = 1:n
    RP_CHECK(i) = Risk_Parity_X(1:n)'* R{i} * Risk_Parity_X(1:n) - Risk_Parity_X(n+1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Market_Weight_r = Market_Weights * mu';
Market_Weight_VAR = Market_Weights * Q * Market_Weights';
Market_Weight_SR = (Market_Weight_r - r_f)/((Market_Weight_VAR)^0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculating Actual Returns in July

July_Returns = zeros(1, Assets_count);
for i = 1:Assets_count
    July_Returns(i) = (Assets_value{i}.AdjClose(903) / Assets_value{i}.AdjClose(883))-1;
end

July_rf = 0.02/12;

MVO_July_r = MVO_X(1,:)* July_Returns';
MVO_July_SR = (MVO_July_r - July_rf)/((MVO_VAR)^0.5);

MVO_R_90_July_r = MVO_R_90_X(1:n)* July_Returns';
MVO_R_90_July_SR = (MVO_R_90_July_r - July_rf)/((MVO_R_90_VAR)^0.5);

MVO_R_95_July_r = MVO_R_95_X(1:n)* July_Returns';
MVO_R_95_July_SR = (MVO_R_95_July_r - July_rf)/((MVO_R_95_VAR)^0.5);

Risk_Parity_July_r = Risk_Parity_X(1:n)'* July_Returns';
Risk_Parity_July_SR = (Risk_Parity_July_r - July_rf)/((Risk_Parity_VAR)^0.5);

Market_Weight_July_r = Market_Weights * July_Returns';
Market_Weight_July_SR = (Market_Weight_July_r - July_rf)/((Market_Weight_VAR)^0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculating Actual Returns in August

August_Returns = zeros(1, Assets_count);
for i = 1:Assets_count
    August_Returns(i) = (Assets_value{i}.AdjClose(926) / Assets_value{i}.AdjClose(903))-1;
end

August_rf = 0.01/12;

MVO_August_r = MVO_X(1,:)* August_Returns';
MVO_August_SR = (MVO_August_r - August_rf)/((MVO_VAR)^0.5);

MVO_R_90_August_r = MVO_R_90_X(1:n)* August_Returns';
MVO_R_90_August_SR = (MVO_R_90_August_r - August_rf)/((MVO_R_90_VAR)^0.5);

MVO_R_95_August_r = MVO_R_95_X(1:n)* August_Returns';
MVO_R_95_August_SR = (MVO_R_95_August_r - August_rf)/((MVO_R_95_VAR)^0.5);

Risk_Parity_August_r = Risk_Parity_X(1:n)'* August_Returns';
Risk_Parity_August_SR = (Risk_Parity_August_r - August_rf)/((Risk_Parity_VAR)^0.5);

Market_Weight_August_r = Market_Weights * August_Returns';
Market_Weight_August_SR = (Market_Weight_August_r - August_rf)/((Market_Weight_VAR)^0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculating Actual Returns in September

September_Returns = zeros(1, Assets_count);
for i = 1:Assets_count
    September_Returns(i) = (Assets_value{i}.AdjClose(947) / Assets_value{i}.AdjClose(926))-1;
end

September_rf = 0.02/12;

MVO_September_r = MVO_X(1,:)* September_Returns';
MVO_September_SR = (MVO_September_r - September_rf)/((MVO_VAR)^0.5);

MVO_R_90_September_r = MVO_R_90_X(1:n)* September_Returns';
MVO_R_90_September_SR = (MVO_R_90_September_r - September_rf)/((MVO_R_90_VAR)^0.5);

MVO_R_95_September_r = MVO_R_95_X(1:n)* September_Returns';
MVO_R_95_September_SR = (MVO_R_95_September_r - September_rf)/((MVO_R_95_VAR)^0.5);

Risk_Parity_September_r = Risk_Parity_X(1:n)'* September_Returns';
Risk_Parity_September_SR = (Risk_Parity_September_r - September_rf)/((Risk_Parity_VAR)^0.5);

Market_Weight_September_r = Market_Weights * September_Returns';
Market_Weight_September_SR = (Market_Weight_September_r - September_rf)/((Market_Weight_VAR)^0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculating Actual Returns in October

October_Returns = zeros(1, Assets_count);
for i = 1:Assets_count
    October_Returns(i) = (Assets_value{i}.AdjClose(968) / Assets_value{i}.AdjClose(947))-1;
end

October_rf = 0.01/12;

MVO_October_r = MVO_X(1,:)* October_Returns';
MVO_October_SR = (MVO_October_r - October_rf)/((MVO_VAR)^0.5);

MVO_R_90_October_r = MVO_R_90_X(1:n)* October_Returns';
MVO_R_90_October_SR = (MVO_R_90_October_r - October_rf)/((MVO_R_90_VAR)^0.5);

MVO_R_95_October_r = MVO_R_95_X(1:n)* October_Returns';
MVO_R_95_October_SR = (MVO_R_95_October_r - October_rf)/((MVO_R_95_VAR)^0.5);

Risk_Parity_October_r = Risk_Parity_X(1:n)'* October_Returns';
Risk_Parity_October_SR = (Risk_Parity_October_r - October_rf)/((Risk_Parity_VAR)^0.5);

Market_Weight_October_r = Market_Weights * October_Returns';
Market_Weight_October_SR = (Market_Weight_October_r - October_rf)/((Market_Weight_VAR)^0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculating Actual Returns in November

November_Returns = zeros(1, Assets_count);
for i = 1:Assets_count
    November_Returns(i) = (Assets_value{i}.AdjClose(988) / Assets_value{i}.AdjClose(968))-1;
end

November_rf = 0.01/12;

MVO_November_r = MVO_X(1,:)* November_Returns';
MVO_November_SR = (MVO_November_r - November_rf)/((MVO_VAR)^0.5);

MVO_R_90_November_r = MVO_R_90_X(1:n)* November_Returns';
MVO_R_90_November_SR = (MVO_R_90_November_r - November_rf)/((MVO_R_90_VAR)^0.5);

MVO_R_95_November_r = MVO_R_95_X(1:n)* November_Returns';
MVO_R_95_November_SR = (MVO_R_95_November_r - November_rf)/((MVO_R_95_VAR)^0.5);

Risk_Parity_November_r = Risk_Parity_X(1:n)'* November_Returns';
Risk_Parity_November_SR = (Risk_Parity_November_r - November_rf)/((Risk_Parity_VAR)^0.5);

Market_Weight_November_r = Market_Weights * November_Returns';
Market_Weight_November_SR = (Market_Weight_November_r - November_rf)/((Market_Weight_VAR)^0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function value = objective(vars, R, n)
    x = vars(1:n); % Extracting x from vars
    xn1 = vars(n+1); % Extracting xn+1 from vars
    
    value = 0;
    for i = 1:n
        value = value + (x' * R{i} * x - xn1)^2;
    end
end