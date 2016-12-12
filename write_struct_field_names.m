% function to write struct field names
% By John Rangel 11/29/16

function data = write_struct_field_names(field_names)
    
    % get number of field names:
    number_of_fields = length(field_names);
    
    % create struct fields:
    for i = 1:number_of_fields
        my_struct.(field_names{i}) = [];
    end
    
    % pass result:
    data = my_struct;