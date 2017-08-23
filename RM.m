classdef RM < handle
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%                               LICENSE                             %%% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Copyright 2017 Mario Senden
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published 
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
% 
% You should have received a copy of the GNU Lesser General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%                             DESCRIPTION                           %%% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% rm = RM() creates an instance of the recurrent model of perceptual
% learning in early visual cortex using standard parameter values 
% (see Lange, Senden, Radermacher, De Weerd, submitted).
% 
% The public properties 'Phi' (reference orientation), 'trials' (number 
% of trials per session) and 'OD' (orientation difference)  control the 
% experimental setup while the property 'mean_JND' (average just noticeable
% difference) stores simulation results.
% 
% Use rm.session() to simulate a single session of staircase experiment. 
% Use rm.reset_weights() to restore the weight matrix to its naive state.
% Use rm.fix(P) to fix a proportion 'P' of connection weights.

    properties (Access = public)
        Phi                     % stimulus orientation
        trials                  % number of staircase trials
        OD                      % OD of current trial
        mean_JND                % JND (average over trials)
    end
    properties (Access = private)
        % functions
        dV                      % neuron dynamics       
        Adiff                   % angluar difference (180° range)
        Cprob                   % connection probability
        
        % parameters
        N                       % number of neurons
        sigma_ff                % width of feedforward orientation bias
        alpha                   % gain of spike encoder
        J_ff                    % feedforward connection strength
        J_rec                   % recurrent connection strength
        a_e                     % exponent excitatory connections
        a_i                     % exponent inhibitory connections
        c_e                     % normalization excitatory connections
        c_i                     % normalization inhibitory connections
        k                       % scaling of variance
        C                       % decision criterion
        eta                     % learning rate
        t_sim                   % simulation time
        tau                     % membrane time constant
        Theta                   % preferred orientation of each neuron
        V                       % membrane potential
        OD_0                    % baseline OD
        W_0                     % baseline connectivity    
        W                       % connection weight matrix
        S                       % selected connections
        
        % auxiliary
        counter                 % keeps track of correct responses
        
    end
    
    methods (Access = public)
        % constructor
        function self = RM(varargin) 
            p = inputParser;
            addOptional(p,'N',512);
            addOptional(p,'alpha',10);
            addOptional(p,'sigma_ff',45);
            addOptional(p,'J_ff',.5);
            addOptional(p,'J_rec',1);
            addOptional(p,'a_e',2.2);
            addOptional(p,'a_i',1.4);
            addOptional(p,'c_e',1.2025e-3);
            addOptional(p,'c_i',1.6875e-3);
            addOptional(p,'k',1.47);
            addOptional(p,'C',.53);
            addOptional(p,'eta',1.5e-9);
            addOptional(p,'t_sim',.5);
            addOptional(p,'tau',15e-3);
            addOptional(p,'OD',7.5);
            
            p.parse(varargin{:});
            
            self.dV         = @(V_rec,V_ff,W,alpha,tau)...
                                (-V_rec+V_ff+W*alpha*max(V_rec,0))/tau;
            self.Adiff      = @(A,B)...
                                angle(exp((A-B)*pi/90*1i))*90/pi;
            self.Cprob      = @(In,a,c)...
                                c*(cosd(2*In)+1).^a;
            
            self.N          = p.Results.N;
            self.alpha      = p.Results.alpha;
            self.sigma_ff   = p.Results.sigma_ff;
            self.J_ff       = p.Results.J_ff;
            self.J_rec      = p.Results.J_rec; 
            self.a_e        = p.Results.a_e;
            self.a_i        = p.Results.a_i;
            self.c_e        = p.Results.c_e;
            self.c_i        = p.Results.c_i; 
            self.k          = p.Results.k;
            self.C          = p.Results.C;
            self.eta        = p.Results.eta;
            self.tau        = p.Results.tau;
            self.t_sim      = p.Results.t_sim;
            
            self.OD_0       = p.Results.OD;
            self.OD         = self.OD_0;
            self.Theta      = linspace(-90,90,self.N)';
            self.Phi        = [];
            self.mean_JND   =  0;
            self.W_0        = (self.Cprob(meshgrid(self.Theta)...
                                -meshgrid(self.Theta)',...
                                self.a_e,self.c_e)...
                                -self.Cprob(meshgrid(self.Theta)...
                                -meshgrid(self.Theta)',...
                                self.a_i,self.c_i));
            self.W          = self.W_0;
            self.S          = ones(self.N);
            self.counter    =   0;
            self.trials     = 480; 
        end
        
        % fixing connections
        function fix(self,p)
           self.S           = double(rand(self.N)>p); 
        end
        
        % resetting weights
        function reset_weights(self)
            self.W          = self.W_0;
        end
        
        % resetting orientation difference
        function set_OD(self,varargin)
            p               = inputParser;
            addOptional(p,'od',[]);
            p.parse(varargin{:});
            od              = p.Results.od;
            if ~isempty(od)
                self.OD         = od;
            else
                self.OD         = self.OD_0;
            end
            
        end
            
        % simulation of training session
        function session(self)
            self.V          = zeros(self.N,1);
            self.mean_JND   = 0;
            self.counter    = 0;
            for t=1:self.trials
               self.trial();
               self.mean_JND    = self.mean_JND+...
                                    (self.OD-self.mean_JND)/t;
            end
        end
    end
    
    methods (Access = private)
        % simulation of individual trial
        function trial(self)
            V_ff            = self.J_ff*exp(...
                                -((self.Adiff(self.Theta,self.Phi)).^2)...
                                /(2*self.sigma_ff^2));
            [~,v]           = ode45(@(t,v)self.dV(v,V_ff,...
                                self.W,self.alpha,self.tau),...
                                [0 self.t_sim],self.V);
            self.V          = v(end,:)';
            M_ref           = self.alpha*max(self.V,0)*self.t_sim;  
            Phi_probe       = self.Phi+((rand()>.5)*2-1)*self.OD;
            V_ff            = self.J_ff*exp(...
                                -((self.Adiff(self.Theta,Phi_probe)).^2)...
                                /(2*self.sigma_ff^2));
            [~,v]           = ode45(@(t,v)self.dV(v,V_ff,...
                                self.W,self.alpha,self.tau),...
                                [0 self.t_sim],self.V);              
            self.V          = v(end,:)';
            M_probe         = self.alpha*max(self.V,0)*self.t_sim;                                
            D               = abs(M_ref-M_probe)./...
                                ((self.k*(M_ref+M_probe)).^.5);
            p               = .5*erfc(-D/(2^.5));
            p(isnan(p))     = .5;
            correct         = mean(p>rand(self.N,1))>=self.C;
            
            if ~correct
                dW              = self.eta*self.S.*(self.V*self.V');
                self.W          = self.W-dW;
                self.OD         = self.OD*1.2;
                self.counter    = 1;
            else
                if self.counter==4
                    self.OD         = self.OD/1.2;
                    self.counter    = 1;
                else
                    self.counter    = self.counter+1;
                end
            end                 
        end
    end
end