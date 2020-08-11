%%
%     Handy Anonymous Functions for MATLAB
%     Copyright (C) 2016  Mehran M. Spitmaan | mehran.m.spitmaan@gmail.com
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

%%
%% cell to mat seperator
cell22mat = @(c) c{:};

%% ###############################
% #### Choose Parenteses
% 
% #### Example(s):
% paren(get(0, 'ScreenSize'), 3:4)  # Returns specific part of the results
% 
% paren(magic(3), 1:2, 2:3)         # Returns specific part of the Matrix
% 
parenn = @(x, varargin) x(varargin{:});
% ###############################

%% ###############################
% #### Choose Curl
% 
% #### Example(s):
% curly({fprintf('This is the first line. \n'),...  # 1st Part (a command)
%     fprintf('This is the second line. \n'),...    # 2nd Part (a command)
%     max(rand(1,10))},...                          # 3rd Part (a command) | We want this result     
%     3)                                            # Indicator
% 
% [A, B] =...                                       # Multiple Outputs
%     curly({'AA','BB'},':')                            # Seperate Cells
% 
curly = @(x, varargin) x{varargin{:}};
% ###############################

%% ###############################
% #### If-Then-Else
% 
% #### Example(s):
% iif(1>0, true,...         # If
%     2>3, false)           # Else If
% 
% iif(2>3, false,...        # If
%     5>7, false,...        # Else If
%     true, 'YES')          # Else
% 
iif = @(varargin)varargin{2*find([varargin{1:2:end}],1,'first')}();
% ###############################

%% ###############################
% #### Map Multiple Functions to a single value (Simple output | Outputs should have same type and size)
% 
% #### Example(s):
% tmpMat = ...              # Output is a Matrix
%     map({rand(1,10)},...  # Value | Shuold be feed as cell {}
%     {@mean,...            # Function(s) | Should be feed as cell {} | 1st Funnction (MATLAB mean function)     
%     @std,...              # 2nd Funnction (MATLAB std function)
%     @(x)(sum(x)/55)})     # 3rd Funnction (User-defined function)
% 
map = @(val, fcns) cellfun(@(f) f(val{:}), fcns);
% ###############################

%% ###############################
% #### Map Multiple Functions to a single value (Cell Output | Output could have diffrent type or size)
% 
% #### Example(s):
% tmpCell = ...             # Output is a Cell
%     mapc({rand(1,10)},... # Value | Shuold be feed as cell {}
%     {@mean,...            # Function(s) | Should be feed as cell {} | 1st Funnction (MATLAB mean function)     
%     @std,...              # 2nd Funnction (MATLAB std function)
%     @(x)(sum(x)/55)})     # 3rd Funnction (User-defined function)
% 
mapc = @(val, fcns) cellfun(@(f) f(val{:}), fcns, 'UniformOutput', false);
% ###############################

%% ###############################
% #### Map Multiple Functions to a single value (Seperate Matrices Outputs | Output could
% #### have diffrent type or size)
% 
% #### Example(s):
% [A B C] = ...             # Multiple Outputs
%     mapm({rand(1,10)},... # Value | Shuold be feed as cell {}
%     {@mean,...            # Function(s) | Should be feed as cell {} | 1st Funnction (MATLAB mean function)     
%     @std,...              # 2nd Funnction (MATLAB std function)
%     @(x)(sum(x)/55)})     # 3rd Funnction (User-defined function)
% 
mapm = @(val, fcns) curly(cellfun(@(f) f(val{:}), fcns, 'UniformOutput', false),':');
% ###############################





%% ###############################
% #### Map Multiple Functions to a Multiple values (Cell Output | Output could have diffrent type or size)
% 
% #### Example(s):
% tmpCell = ...             # Output is a Cell
%     mapfc({rand(1,10),... # 1st Value | Shuold be feed as cell {}
%     rand(20,1)},...       # 2nd Value
%     {@mean,...            # Function(s) | Should be feed as cell {} | 1st Funnction (MATLAB mean function)     
%     @std,...              # 2nd Funnction (MATLAB std function)
%     @(x)(sum(x)/55)})     # 3rd Funnction (User-defined function)
%
mapfc = @(val, fcns) cellfun(@(f) cellfun(f,val,'UniformOutput',false), fcns,'UniformOutput',false);
% ###############################

%% ###############################
% #### Map Multiple Functions to a Multiple values (Seperate Matrices Outputs | Output could
% #### have diffrent type or size)
% 
% #### Example(s):
% [A B C] = ...             # Multiple Outputs
%     mapfm({rand(1,10),... # 1st Value | Shuold be feed as cell {}
%     rand(20,1)},...       # 2nd Value
%     {@mean,...            # Function(s) | Should be feed as cell {} | 1st Funnction (MATLAB mean function)     
%     @std,...              # 2nd Funnction (MATLAB std function)
%     @(x)(sum(x)/55)})     # 3rd Funnction (User-defined function)
%
mapfm = @(val, fcns) curly(cellfun(@(f) cellfun(f,val,'UniformOutput',false), fcns,'UniformOutput',false),':');
% ###############################





%% ###############################
% #### Recursion Function
% 
% #### Example(s):
% factorial = @(n)...
%     recur(@(f,k)...
%         iif(k == 0, 1,...
%         true, @() k * f(f,k-1)),...
%     n)
% factorial(5)
% 
recur = @(f, varargin) f(f, varargin{:});
% ###############################

%% ###############################
% #### Loop Iterator
% 
% #### Example(s):
% factorial = @(n) loop({1, 1}, ...     # Start with k = 1 and x = 1
%     @(k,x) k <= n, ...                # While k <= n
%     @(k,x) {k + 1, ...                #   k = k + 1;
%     k * x}, ...                       #   x = k * x;
%     @(k,x) x);                        # End, returning x
% factorial(5)
% 
loop  = @(x0,c,f,r)recur(@(g,x)iif(c(x{:}),@()g(g,f(x{:})),1,r(x{:})),x0);
% ###############################





%% ############################### 
% #### Pairwise function-value (Cell Output | Output could have diffrent type or size)
% 
% #### Example(s):
% tmpCell = ...                         # Output is a Cell
% mappc({@max, rand(1,20),...            # 1st Function-Value Pair | Shuold be feed as cell {}
%     @mean, rand(40,1),...             # 2nd Function-Value Pair
%     @(x)sum(x)/55, rand(10)})         # 3rd Function-Value Pair
% 
mappc = @(v) loop({1,{}},@(a,b) a <= size(v,2),@(a,b) {a+2,{b{:} v{a}(v{a+1})}},@(a,b) b);
% ###############################

%% ############################### 
% #### Pairwise function-value (Seperate Output | Output could have diffrent type or size)
% 
% #### Example(s):
% [A B C] = ...             # Multiple Outputs
% mappm({@max, rand(1,20),...            # 1st Function-Value Pair | Shuold be feed as cell {}
%     @mean, rand(40,1),...             # 2nd Function-Value Pair
%     @(x)sum(x)/55, rand(10)})         # 3rd Function-Value Pair
% 
mappm = @(v) curly(loop({1,{}},@(a,b) a <= size(v,2),@(a,b) {a+2,{b{:} v{a}(v{a+1})}},@(a,b) b),':');
% ###############################

%% ###############################
% #### Matrix Initializer 
% 
% #### Example(s):
% [A1,A2,A3] = ...          # Multiple Outputs
%     matIni(3);            # Calling Initializer
% 
matIni = @(N) cell22mat(cellfun(@(x) [],mat2cell(zeros(1,N),1,ones(1,N)),'UniformOutput',false)); 
% ###############################


%% warning('Anonymous Functions Loaded!!');