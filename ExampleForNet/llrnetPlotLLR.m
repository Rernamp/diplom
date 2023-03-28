function varargout = llrnetPlotLLR(simResults,layoutTitle)
%llrnetPlotLLR Plot LLR curves for QAM symbols
%   llrnetPlotLLR(R,TITLE) plots LLR versus real part of the received
%   symbol for QAM symbols based on the results structure, R. This function
%   plots the LLR values for the odd bits. TITLE is the title of the plot
%   and must be a string or character array.
%
%   R is generated by the llrnetQAMLLR function and has the following
%   fields:
%
%   exactLLR        - Exact LLR
%   approxLLR       - Approximate LLR
%   predictedLLR    - LLRNet predicted LLR
%   RxSymbols       - Received symbols
%   M               - Modulation order
%   SNRValues       - SNR values as a scalar or row vector
%   HiddenLayerSize - Hidden layer size for neural network
%   NumSymbols      - Number of symbols to simulate
%
%   See also llrnetQAMLLR, LLRNeuralNetworkExample.

%   Copyright 2019 The MathWorks, Inc.

exactLLR = simResults.exactLLR;
approxLLR = simResults.approxLLR;
predictedLLR = simResults.predictedLLR;
M = simResults.M;
SNRValues = simResults.SNRValues;
rxSym = simResults.RxSymbols;

xLimitMax = ceil(max(real(rxSym(:))));
xLimitMin = floor(min(real(rxSym(:))));

k = log2(M);

fig = figure;
set(fig, 'Position', [375 300 990 745]);
t = tiledlayout(k/2,length(SNRValues));
for p=1:2:k
  
  for snrIdx = 1:length(SNRValues)
    if mod(p,2) == 1
      symbols = real(rxSym(:,snrIdx));
      xLabelText = 'Re(Symbol)';
    else
      symbols = imag(rxSym(:,snrIdx));
      xLabelText = 'Im(Symbol)';
    end
    
    indices = find(abs(symbols) <= max(abs([xLimitMin xLimitMax])));
    
    exactLLRSNR = reshape(exactLLR(:,snrIdx), k, []);
    maxLogLLRSNR = reshape(approxLLR(:,snrIdx), k, []);
    predictedLLRSNR = reshape(predictedLLR(:,snrIdx), k, []);
    
    exactLLRSNR = exactLLRSNR(:,indices);
    maxLogLLRSNR = maxLogLLRSNR(:,indices);
    predictedLLRSNR = predictedLLRSNR(:,indices);
    
    yLimitMax = max([exactLLRSNR(p:k:end,:) maxLogLLRSNR(p:k:end,:) predictedLLRSNR(p:k:end,:)]);
    yLimitMin = min([exactLLRSNR(p:k:end,:) maxLogLLRSNR(p:k:end,:) predictedLLRSNR(p:k:end,:)]);
    
    nexttile
    plot(symbols(indices,1), predictedLLRSNR(p:k:end,:), 'o')
    hold on
    plot(symbols(indices,1), exactLLRSNR(p:k:end,:), '.')
    plot(symbols(indices,1), maxLogLLRSNR(p:k:end,:), '.')
    hold off
    axis([xLimitMin xLimitMax yLimitMin yLimitMax])
    grid on
    xlabel(xLabelText)
    ylabel(sprintf('LLR of bit #%d',p))
    if p==1
      title(sprintf('SNR=%ddB',SNRValues(snrIdx)))
      if snrIdx == 1
        legend('LLRNet','Exact LLR','Approx. LLR','Location','northwest')
      end
    end
  end
end
title(t, layoutTitle)

if nargout > 0
  varargout{1} = fig;
end
end
