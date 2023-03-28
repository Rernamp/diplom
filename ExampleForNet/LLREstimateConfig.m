classdef LLREstimateConfig
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        amplifaer
        useNet
    end

    methods
        function obj = LLREstimateConfig(useNet, amplifaer)
            obj.amplifaer = amplifaer;
            obj.useNet = useNet;
        end
        
        function need = UseNet(obj)
            need = obj.useNet;
        end
        
        function outputSignal = amplifaerSignal(obj, input)
            outputSignal = obj.amplifaer(input);
        end
    end
end