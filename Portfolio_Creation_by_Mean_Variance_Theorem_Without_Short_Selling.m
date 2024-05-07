%clear
clearvars -except std_devi_without_3 goal_R_without_3
clc

%define the series of CSV file names
ETF = {'SPY.csv', 'GOVT.csv', 'EEMV.csv', 'CME.csv', 'BR.csv', 'CBOE.csv', 'ICE.csv', 'ACN.csv'};
ETF_count = numel(ETF);

%define a cell array to store the data from each CSV file
ETF_value = cell(1, ETF_count);
%loop through each CSV file
for i = 1:ETF_count
    %read the CSV file
    ETF_value{i} = readtable(ETF{i});
    %convert the first column to datetime format
    ETF_value{i}.Date = datetime(ETF_value{i}.Date);
    %extract month and year information
    ETF_value{i}.Month = month(ETF_value{i}.Date);
end

%define a cell array to store return of each end of month
end_of_month = cell(1, ETF_count);
%getting total number of months
monthYear = year(ETF_value{1}.Date) * 100 + month(ETF_value{1}.Date);
total_months_count = size(histcounts(monthYear, unique(monthYear)),2);
zeros_1 = zeros(total_months_count,1);

%getting the end of month values
for i = 1:ETF_count
    past_month = 12;
    p = 1;
    end_of_month{i} = table(zeros_1, zeros_1, 'VariableNames', {'End_of_Month', 'Return'});
    for j = 1:height(ETF_value{i})
        current_month = ETF_value{i}.Month(j);
        if current_month ~= past_month
            end_of_month{i}.End_of_Month(p) = past_month;
            end_of_month{i}.Return(p) = ETF_value{i}.AdjClose(j-1);
            p = p+1;
        end
        past_month = ETF_value{i}.Month(j);
    end
end

%define a cell array to store months returns
months_returns = cell(1, ETF_count);
%define a matrix to store arithmatic average
Average_returns = zeros(1, ETF_count);
%define a matrix to store ETF returns
ETF_returns = zeros(1, ETF_count);
%getting included months number (T)
months_count = total_months_count - 1;
zeros_2 = zeros(months_count,1);


for i = 1:ETF_count
    months_returns{i} = table(zeros_2, zeros_2, zeros_2, 'VariableNames', {'Month', 'Return', 'Return_minus_Average'});
    for j = 2:height(end_of_month{i})
        months_returns{i}.Month(j-1) = end_of_month{i}.End_of_Month(j);
        %calculating monthly_returns(i) = (last_day_of_month(i)/ last_day_of_month(i-1))-1
        months_returns{i}.Return(j-1) = ((end_of_month{i}.Return(j)/end_of_month{i}.Return(j-1))-1);
    end
    %calcualting arithmatic average = avg(r_it)
    Average_returns(i) = mean(months_returns{i}.Return);
    %calcualting geometric expected return of asset i = ((product(1+r_it))^(1/T))-1
    ETF_returns(i) = ((prod((months_returns{i}.Return)+1))^(1/months_count))-1;
    
    %calculating r_it - avg(r_it)
    for k = 1:months_count
        months_returns{i}.Return_minus_Average(k) = months_returns{i}.Return(k)-Average_returns(i);
    end
end

%define a matrix to store ETF covariance
ETF_covariance = zeros(ETF_count, ETF_count);
for i = 1:ETF_count
    for j = 1:ETF_count
        Difference_product = 1;
        Difference_product_sum = 0;
        if ETF_covariance(i,j) == 0 %to avoid calculating the covariance twice
            for k = 1:months_count
                %[r_it - avg(r_it)] * [r_jt - avg(r_jt)]
                Difference_product = months_returns{i}.Return_minus_Average(k) * months_returns{j}.Return_minus_Average(k);
                %sum{[r_it - avg(r_it)] * [r_jt - avg(r_jt)]}
                Difference_product_sum = Difference_product_sum + Difference_product;
            end
            %covariance between assets i and j = sum{[r_it - avg(r_it)] * [r_jt - avg(r_jt)]}/T
            ETF_covariance(i,j) = Difference_product_sum/months_count;
            ETF_covariance(j,i) = ETF_covariance(i,j);
        end
    end
end



n=ETF_count;
mu=ETF_returns; % expected returns of assets

Q=ETF_covariance; % covariance matrix
c=zeros(n,1); %linear coefficients
A = -mu; %unequal A matrix
Aeq =[ones(1,n)]; %equal A matrix
beq =[1]; %equal b vector
ub = [inf; inf; inf; inf; inf; inf; inf; inf;]; %upper bound
lb_without = [0; 0; 0; 0; 0; 0; 0; 0;]; %lower bound without short selling
lb_5 = [0.05; 0.05; 0.05; 0.05; 0.05; 0.05; 0.05; 0.05;]; %lower bound 5%

%compute minimum variance portfolio:
%without shorting
[x_min_without(1,:),fval_min_without(1,1)] = quadprog (Q, c, [], [], Aeq, beq, lb_without, ub);
r_min_without = x_min_without(1,:)*mu';
%at least 5%
[x_min_5(1,:),fval_min_5(1,1)] = quadprog (Q, c, [], [], Aeq, beq, lb_5, ub);
r_min_5 = x_min_5(1,:)*mu';

%compute maximum variance portfolio for "at least 5%"
f = -ETF_returns;
[x_max_5(1,:),fval_max_5(1,1)] = linprog(f, [], [], Aeq, beq, lb_5, ub);

%expected return goals range from minimum variance portfolio to maximum return between assets
%without shorting
goal_R_without = linspace(r_min_without, max(ETF_returns), 10);
%at least 5%
goal_R_5 = linspace(r_min_5, -fval_max_5(1,1), 10);

%efficient frontier values without shorting
for a=1:length(goal_R_without)
    b = -goal_R_without(a);
    [x_without(a,:),fval_without(a,1)] = quadprog (Q, c, A, b, Aeq, beq, lb_without, ub);
    std_devi_without(a,1)=(fval_without(a,1)*2)^0.5; %calculating standard deviation from objective function value
end

%efficient frontier values at least 5%
for a=1:length(goal_R_5)
    b = -goal_R_5(a);
    [x_5(a,:),fval_5(a,1)] = quadprog (Q, c, A, b, Aeq, beq, lb_5, ub);
    std_devi_5(a,1)=(fval_5(a,1)*2)^0.5; %calculating standard deviation from objective function value
end


hold on
plot(std_devi_without, goal_R_without, '-k*');
plot(std_devi_5, goal_R_5, '-ro');
xlabel('volatility sigma')
ylabel('expected return goal R')
title('The efficient frontier of MVO - 8 Assets')
grid on;
legend('Without\_Shorting', 'At least\_5%', 'Location', 'best');
hold off

hold on
plot(std_devi_without_3, goal_R_without_3, '-b^');
plot(std_devi_without, goal_R_without, '-k*');
plot(std_devi_5, goal_R_5, '-ro');
xlabel('volatility sigma')
ylabel('expected return goal R')
title('The efficient frontier of MVO - Mixed')
grid on;
legend('Without\_Shorting\_3 Assets','Without\_Shorting\_8 Assets', 'At least\_5%\_8 Assets', 'Location', 'best');
hold off