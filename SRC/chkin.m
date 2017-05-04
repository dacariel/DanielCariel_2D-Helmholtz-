function s = chkin(s, field)
% Purpose:
% Adds or creates one or more empty fields specified by cell array "field"
% Example: >> chkin(s, {'a'});                % check in 1 item 

for k=1:length(field)
   s.(field{k}) = [];
end

end   