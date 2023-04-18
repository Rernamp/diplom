close all;
subsystemType = '32APSK 4/5'; %#ok<UNRCH>
EsNoValues = 13.2:0.1:13.6;     % in dB
numFrames = 1;
numErrors = 200;

trainNow = true;
if trainNow && (~strcmp(subsystemType,'16APSK 2/3') || ~isequal(EsNoValues,8.6:0.1:9))
    % Train the networks for each EsNo value
    numTrainSymbols = 1e2;
    hiddenLayerSize = 8;
    llrNets = llrnetTrainDVBS2LLRNetwork(subsystemType, EsNoValues, numTrainSymbols, hiddenLayerSize);
else
    load('llrnetDVBS2Networks','llrNets','subsystemType','EsNoValues');
end

estimateConfig = LLREstimateConfig(true, @(input) (amplifaerSaleh(input)));

% Simulate PER with exact LLR, approximate LLR, and LLRNet
[perLLRSalef,perApproxLLRSalef,perLLRNetSalef] = customLlrnetDVBS2PER(subsystemType,EsNoValues,llrNets,numFrames,numErrors, estimateConfig);
type = subsystemType + " + Saleh";
llrnetPlotLLRvsEsNo(perLLRSalef,perApproxLLRSalef,perLLRNetSalef, EsNoValues,type)

%%estimateConfig = LLREstimateConfig(false, @(input) (amplifaerGhorbani(input)));

% Simulate PER with exact LLR, approximate LLR, and LLRNet
%[perLLRGhorbani,perApproxLLRGhorbani,perLLRNetGhorbani] = customLlrnetDVBS2PER(subsystemType,EsNoValues,llrNets,numFrames,numErrors, estimateConfig);
%type = subsystemType + " + Ghorbani";
%llrnetPlotLLRvsEsNo(perLLRGhorbani,perApproxLLRGhorbani,perLLRNetGhorbani,EsNoValues,type)

estimateConfig = LLREstimateConfig(true, @(input) (amplifaerDummy(input)));

[perLLRDummy,perApproxLLRDummy,perLLRNetDummy] = customLlrnetDVBS2PER(subsystemType,EsNoValues,llrNets,numFrames,numErrors, estimateConfig);
type = subsystemType;
llrnetPlotLLRvsEsNo(perLLRDummy,perApproxLLRDummy,perLLRNetDummy,EsNoValues,type)

save("result.mat","perLLRNetDummy","perLLRDummy","perApproxLLRDummy", "perLLRNetSalef","perApproxLLRSalef","perLLRSalef");

figure ()
semilogy(EsNoValues, perLLRDummy(:,1))
hold on
semilogy(EsNoValues, perLLRSalef(:,1))
xlabel('E_s/N_o (dB)')
ylabel('PER')
grid on
legend('Exact LLR', 'Exact Saleh LLR')

figure ()
semilogy(EsNoValues, perLLRDummy(:,1))
hold on
semilogy(EsNoValues, perLLRSalef(:,1))
semilogy(EsNoValues, perLLRNetDummy(:,1))
semilogy(EsNoValues, perLLRNetSalef(:,1))
xlabel('E_s/N_o (dB)')
ylabel('PER')
grid on
legend('Exact LLR', 'Exact Saleh LLR', 'LLR Net', 'LLR Net Saleh')
% figure ()
% semilogy(EsNoValues, perLLRDummy(:,1))
% hold on
% semilogy(EsNoValues, perLLRGhorbani(:,1))
% semilogy(EsNoValues, perApproxLLRDummy(:,1))
% semilogy(EsNoValues, perApproxLLRGhorbani(:,1))
% xlabel('E_s/N_o (dB)')
% ylabel('PER')
% grid on
% legend('Exact LLR', 'Exact Ghorbani LLR','Approx. LLR', 'Approx. Ghorbani LLR')
% 
% function output = amplifaerGhorbani(input)

% x_1 = 1.92;
% x_2 = 1.74;
% x_3 = 0.92;
% x_4 = 0;
% x_5 = 1.74;
% y_1 = 0.02;
% y_2 = 1;
% y_3 = 0.4;
% y_4 = 0;
% y_5 = 3.5;
% 
% G = (x_1.*abs(xComplex).^x_2)./(1+x_3.*abs(xComplex).^x_5)+x_4.*abs(xComplex);
% F = (y_1.*abs(xComplex).^y_2)./(1+y_3.*abs(xComplex).^y_5)+y_4.*abs(xComplex).^y_2;
% 
% output = G .* exp(1i * 2 * pi .* F);
% output = input .* (output ./ abs(input));

% end


function output = amplifaerSaleh(input)

a_A = 2.1587; 
b_A = 1.1517;
a_F = 4.0033;
b_F = 9.1041;

input = input .* 2;

G = a_A.*abs(input)./(1+b_A.*(abs(input).^2));
F = a_F*(abs(input).^2)./(1+b_F*(abs(input).^2));
output = G .* exp(1i * 2 * pi .* F);
output = input .* (output ./ abs(input));

end

function output = amplifaerDummy(input)
    output = input;
end