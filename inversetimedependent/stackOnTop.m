%% function stackOnTop

%
% takes a field name of the structure array and stacks all those 
% vectors on top of each other
%

function output = stackOnTop( struct, fieldname )

    
        aux = arrayfun( @(s) reshape(s.(fieldname), length(s.(fieldname)), 1), struct, ...
            'UniformOutput', 0 );
        output = cat(1, aux{:} );

end