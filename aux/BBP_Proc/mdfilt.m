classdef mdfilt
  methods
% this file contains the median filter functions from Dal'Ollmio et al, 2022
% function for One-dimensional median filter
    function [filter_data] =  mdfilt1(obj,data, kernel_size);
      halfkernel = floor(kernel_size/2);

      filtered_data = NaN.*ones(size(data));

      for n = 1:(length(data));
        i1 = nanmax([1, n-halfkernel]);
        i2 = nanmin([length(data), n+halfkernel]);
        filter_data(n) = nanmedian(data(i1:i2));
      end
    end %function
% function to define adaptive median filtering based on Christina Schallemberg's suggestion for CHLA
    function [ymf] = adaptive_mdfilt1(obj,x, y,win);
    % compute x resolution
      xres = diff(x);

    % initialise medfiltered array
      ymf = ones(size(y)).*NaN;

      ir_GT0 = find(xres>0);
      if any(ir_GT0);
        win_GT0 = win;
        ymf(ir_GT0) = obj.mdfilt1(y(ir_GT0), win_GT0);


      end
    end %function
  end % methods
end %classdef
