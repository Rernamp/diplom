classdef LLREstimateConfig
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        amplifaer
        useNet
        useApproximate
    end

    methods
        function obj = LLREstimateConfig(useNetAndApproximate, amplifaer)
            obj.amplifaer = amplifaer;
            obj.useNet = useNetAndApproximate;
            obj.useApproximate = useNetAndApproximate;
        end
        
        function need = UseNet(obj)
            need = obj.useNet;
        end

        function need = UseApporximate(obj)
            need = obj.useApproximate;
        end
        
        function outputSignal = amplifaerSignal(obj, input)
            outputSignal = obj.amplifaer(input);
        end
    end
end