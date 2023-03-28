
subsystemType = '16APSK 2/3'; %#ok<UNRCH>
EsNoValues = 8.6:0.1:8.9;     % in dB
numFrames = 1000;
numErrors = 200;

trainNow = false;
if trainNow && (~strcmp(subsystemType,'16APSK 2/3') || ~isequal(EsNoValues,8.6:0.1:9))
    % Train the networks for each EsNo value
    numTrainSymbols = 1e4;
    hiddenLayerSize = 64;
    llrNets = llrnetTrainDVBS2LLRNetwork(subsystemType, EsNoValues, numTrainSymbols, hiddenLayerSize);
else
    load('llrnetDVBS2Networks','llrNets','subsystemType','EsNoValues');
end

estimateConfig = LLREstimateConfig(false, @(input) (amplifaer(input)));

% Simulate PER with exact LLR, approximate LLR, and LLRNet
[perLLR,perApproxLLR,perLLRNet] = customLlrnetDVBS2PER(subsystemType,EsNoValues,llrNets,numFrames,numErrors, estimateConfig);
llrnetPlotLLRvsEsNo(perLLR,perApproxLLR,perLLRNet,EsNoValues,subsystemType)


function output = amplifaer(input)


a_A = 2.1587; 
b_A = 1.1517;
a_F = 4.0033;
b_F = 9.1041;

G = a_A.*abs(input)./(1+b_A.*(abs(input).^2));
F = a_F*(abs(input).^2)./(1+b_F*(abs(input).^2));
output = G .* exp(1i * 2 * pi .* F);
output = input .* (output ./ abs(input));

end