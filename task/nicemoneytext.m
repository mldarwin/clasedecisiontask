function formattedtext = nicemoneytext(dollars)
% Function to nicely format money quantities for experimental presentation.
% Formats losses as -$X and gains or zero as $X.
%
% PSH 8/9/21

if nargin < 1
    error('Must provide a dollar quantity.')
elseif nargin > 1
    error('Too many input arguments.')
end
    
if ~isfinite(dollars)
    error('Not a finite quantity of dollars!')
elseif ~isreal(dollars)
    error('Not a real quantity of dollars!')
end


if dollars < 0
    formattedtext = sprintf('-$%0.2f',abs(dollars));
else
    formattedtext = sprintf('$%0.2f',abs(dollars)); % use abs() to account for the edge case in which someone passes "-0" as the dollar quantity.
end