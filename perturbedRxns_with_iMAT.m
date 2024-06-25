clc;
clear all;

initCobraToolbox;
changeCobraSolver('gurobi','all');

model = readCbModel('iAB5075.xml');
model.ub(1246) = 0; % Set upper bound of leakage reaction reaction to 0 
model.lb(1246) = 0; % Set lower bound of of leakage reaction reaction to 0 
model.lb(2140) = 8.39; % Set lower bound of NGAM reaction to 8.39
model.ub(2140) = 8.39; % Set upper bound of NGAM reaction to 8.39
model.lb(2139) = 0.1; % Set lower bound of biomass reaction to 0.1 
model.ub(2026) = -0.01; % Set upper bound of oxygen uptake  reaction to -0.01

[selExc, selUpt] = findExcRxns(model);
uptake = [num2cell(find(selUpt)),model.rxns(selUpt),model.rxnNames(selUpt), num2cell(model.lb(selUpt)),num2cell(model.ub(selUpt))];

[hiCarbonRxns, zeroCarbonRxns, nCarbon] = findCarbonRxns(model, 1); % Find reactions containing at least 1 carbon
idx_zeroC = cell2mat(uptake(ismember(uptake(:,2), zeroCarbonRxns) == 1, 1)); % Index of uptake reactions containing zero carbon metabolites
idx_carbon = cell2mat(uptake(ismember(uptake(:,2), setdiff(model.rxns(selUpt),zeroCarbonRxns)) == 1, 1)); % Index of uptake reactions containing carbon metabolites

model.lb(idx_carbon) = -2.5; % Set lower bounds of carbon-containing uptake reactions to -2.5
model.lb(1958) = -10; % Set lower bound of glucose uptake reaction to -10

model.c(2139) = 1; 
gmodelWT.A = sparse(model.S);
gmodelWT.obj = model.c;
gmodelWT.rhs = model.b;
gmodelWT.lb = model.lb;
gmodelWT.ub = model.ub;
gmodelWT.vtype = 'C';
gmodelWT.sense = '=';
gmodelWT.modelsense = 'max';
fba_iab5075 = gurobi(gmodelWT);
fba_iab5075.x(2139) % Biomass production: 4.09

strains = {'1207552', '1428368', '1457504', '34654', '478810'};

for s = 1:length(strains)
    strain = strains{s};
    
    counts = readtable([strain '_getmm_cpmThreshold_0_5.xlsx']);
    ab_ReadCounts = table2array(counts(:, [2:end-2])); % Exclude last two columns (NACL)

    expressionData.gene = counts.Var1;
    expressionRxns = [];

    % Map expression data to reactions
    for i = 1:size(ab_ReadCounts,2)
        expressionData.value = ab_ReadCounts(:,i);
        [expressionRxns(:,i), ~, ~] = mapExpressionToReactions(model, expressionData);
    end

    % Calculate quartiles for later use
    Q1 = quantile(mean(expressionRxns,2), 0.25);
    Q3 = quantile(mean(expressionRxns,2), 0.75);

    met = [];
    rxn = [];
    models = [];
    genes = [];

    for i = 1:size(ab_ReadCounts,2)
        expressionData.value = ab_ReadCounts(:,i);
        [expressionRxns, ~, gene_used] = mapExpressionToReactions(model, expressionData);
        [tissModel] = iMAT(model, expressionRxns, Q1, Q3);
        met = [met length(tissModel.mets)];
        rxn = [rxn length(tissModel.rxns)];
        lia = ismember(model.rxns,tissModel.rxns);
        models = [models lia];
        genes = [genes gene_used'];
    end

    t = [model.rxns, model.rxnNames, num2cell(models), genes];
    t = cell2table(t);

    % Define and rename columns dynamically
    binaryColumnNames = {'AMI25_1', 'AMI25_2', 'AMI75_1', 'AMI75_2', 'CIP25_1', 'CIP25_2', 'CIP75_1', 'CIP75_2', ...
        'COL25_1', 'COL25_2', 'COL75_1', 'COL75_2', 'MERO25_1', 'MERO25_2', 'MERO75_1', 'MERO75_2', 'MHB_1', 'MHB_2', ...
        'AMI25_1G', 'AMI25_2G', 'AMI75_1G', 'AMI75_2G', 'CIP25_1G', 'CIP25_2G', 'CIP75_1G', 'CIP75_2G', ...
    'COL25_1G', 'COL25_2G', 'COL75_1G', 'COL75_2G', 'MERO25_1G', 'MERO25_2G', 'MERO75_1G', 'MERO75_2G', 'MHB_1G', 'MHB_2G'};
    
    newColumnNames = cellfun(@(x) [strain '_' x], binaryColumnNames, 'UniformOutput', false);
    newColumnNames = [{'Rxns', 'RxnNames'}, newColumnNames];

    t.Properties.VariableNames = newColumnNames;

    % Display and write tables
    disp(t);
    writetable(t, [strain '_iMATModels_binary_GETMM_CPMThreshold_0_5_withUsedGenes.xlsx']);

    % Calculate sums and differences
    binaryColumnNames = t.Properties.VariableNames(3:20);

    sumTable = t(:,1);

    for i = 1:2:length(binaryColumnNames)
        col1 = binaryColumnNames{i};
        col2 = binaryColumnNames{i+1};
        sumValues = sum(t{:, {col1, col2}}, 2);
        newColName = strrep(col1, '_1', '');
        sumTable.(newColName) = sumValues;
    end

    diffTable = table(sumTable.Rxns);
    % Get the column names for the binary columns
    binaryColumnNames = sumTable.Properties.VariableNames(2:end);

    for i = 1:length(binaryColumnNames)
        col1 = binaryColumnNames{i};
        col2 = [strain '_MHB'];
        diffValues = sumTable.(col1) - sumTable.(col2);
        newColName = [col1, ' - ', col2];
        diffTable.(newColName) = diffValues;
    end

    disp(diffTable);
    writetable(diffTable, [strain '_iMATModels_binary_GETMM_CPMThreshold_0_5_diffTable.xlsx']);
end

% Clear unnecessary variables
clear ab_ReadCounts expressionData expressionRxns gene_used genes Q1 Q3 met rxn models newColName newColumnNames sumTable diffTable diffValues counts
