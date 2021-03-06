%==========================================================================
%********************FUNCTION: random_vector(n,deg)************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% n: Size of the vector
% deg: the number of non-zero elements in the vector
% -------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% V: The output random vector
% -------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function generates a random vector with certin length and a given 
% number of 1's(on average).
%--------------------------------------------------------------------------
%==========================================================================
%==========================================================================

function V = random_vector(n,deg)
a=clock;                                % Initialize the seed for random number generation with the clock value.
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*a))); 
%-----------------------------Initialization-------------------------------
V = zeros(1,n);
%--------------------------------------------------------------------------

%---------------------------Generate the Vector----------------------------
for i = 1:n
    p = rand;           % Pick a random number between 0 and 1
    
    %--With probability (deg/n) assign 1 to elements of the output vector--
    if (p <= (deg/n))   
        V(i) = 1;      
    end    
    %----------------------------------------------------------------------
    
end
%--------------------------------------------------------------------------
