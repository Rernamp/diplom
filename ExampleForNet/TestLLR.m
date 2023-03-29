
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

subsystemType = '16APSK 2/3'; %#ok<UNRCH>
EsNoValues = 8.6:0.1:8.9;     % in dB

estimateConfig = LLREstimateConfig(false, @(input) (amplifaerSaleh(input)));

% Simulate PER with exact LLR, approximate LLR, and LLRNet
[perLLR,perApproxLLR,perLLRNet] = customLlrnetDVBS2PER(subsystemType,EsNoValues,llrNets,numFrames,numErrors, estimateConfig);
type = subsystemType + " + Saleh";
llrnetPlotLLRvsEsNo(perLLR,perApproxLLR,perLLRNet,EsNoValues,type)

estimateConfig = LLREstimateConfig(false, @(input) (amplifaerGhorbani(input)));

% Simulate PER with exact LLR, approximate LLR, and LLRNet
[perLLR,perApproxLLR,perLLRNet] = customLlrnetDVBS2PER(subsystemType,EsNoValues,llrNets,numFrames,numErrors, estimateConfig);
type = subsystemType + " + Ghorbani";
llrnetPlotLLRvsEsNo(perLLR,perApproxLLR,perLLRNet,EsNoValues,type)


function output = amplifaerGhorbani(input)

x_1 = 1.92;
x_2 = 1.74;
x_3 = 0.92;
x_4 = 0;
x_5 = 1.74;
y_1 = 0.02;
y_2 = 1;
y_3 = 0.4;
y_4 = 0;
y_5 = 3.5;

G = (x_1.*abs(xComplex).^x_2)./(1+x_3.*abs(xComplex).^x_5)+x_4.*abs(xComplex);
F = (y_1.*abs(xComplex).^y_2)./(1+y_3.*abs(xComplex).^y_5)+y_4.*abs(xComplex).^y_2;

output = G .* exp(1i * 2 * pi .* F);
output = input .* (output ./ abs(input));

end


function output = amplifaerSaleh(input)

a_A = 2.1587; 
b_A = 1.1517;
a_F = 4.0033;
b_F = 9.1041;

G = a_A.*abs(input)./(1+b_A.*(abs(input).^2));
F = a_F*(abs(input).^2)./(1+b_F*(abs(input).^2));
output = G .* exp(1i * 2 * pi .* F);
output = input .* (output ./ abs(input));

end