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
    % see: Lange G, Senden M, Radermacher A, De Weerd P.
    % Interfering with a memory without disrupting its trace (submitted).
    %
    % Use rm.set_OD(x) to set orientation difference to value 'x'; if no value
    %     is provided, OD will be reset to its baseline state (7.5 unless
    %     specified otherwise during construction)
    % Use rm.set_PHI(x) to set reference orientation to value 'x'; if no value
    %     is provided, PHI wil be reset to ts baseline state (135° unless
    %     specified otherwise during construction)
    % Use rm.fix(P) to fix a proportion 'P' of connection weights.
    % Use rm.get_JND() to read out the current JND
    % Use rm.get_weights() to retrieve excitatory and inhibitory weights
    % Use rm.get_response() to read membrane potentials and firing rates
    % Use rm.session() to simulate a single session of staircase experiment.
    % Use rm.reset() to restore the model to its naive state.
    
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
        eta                     % gobal learning rate
        Eta                     % learning rate per connection
        mu                      % exponent of power law weight dependence
        t_sim                   % simulation time
        tau                     % membrane time constant
        Theta                   % preferred orientation of each neuron
        V_0                     % baseline membrane potential
        W_exc                   % excitatory lateral connectivity
        W_inh                   % inhibitory lateral connectivity
        
        
        % experimental setup
        Phi_0                   % baseline stimulus orientation
        Phi                     % current stimulus orientation
        trials                  % number of staircase trials
        OD_0                    % baseline OD
        OD                      % OD of current trial
        p_correct               % probability of receiving "correct"
                                % feedback on incorrect trials
        p_incorrect             % probability of receiving "incorrect"
                                % feedback on correct trials
        mean_JND                % JND (average over trials)
        
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
            addOptional(p,'k',4);
            addOptional(p,'C',.53);
            addOptional(p,'eta',1.5e-11);
            addOptional(p,'mu',0);
            addOptional(p,'t_sim',.5);
            addOptional(p,'tau',15e-3);
            addOptional(p,'trials',480);
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
            self.Eta        = self.eta*ones(self.N);
            self.mu         = p.Results.mu;
            self.tau        = p.Results.tau;
            self.t_sim      = p.Results.t_sim;
            self.trials     = p.Results.trials;
            self.Theta      = linspace(-90,90,self.N)';
            self.mean_JND   =  0;
            self.V_0        = zeros(self.N,1);
            self.W_exc      = self.J_rec * self.Cprob(meshgrid(self.Theta)...
                -meshgrid(self.Theta)',...
                self.a_e,self.c_e);
            self.W_inh      = self.J_rec * self.Cprob(meshgrid(self.Theta)...
                -meshgrid(self.Theta)',...
                self.a_i,self.c_i);
            self.Phi_0      = 135;
            self.Phi        = self.Phi_0;
            self.OD_0       = p.Results.OD;
            self.OD         = self.OD_0;
            self.p_correct  =  0;
            self.p_incorrect=  0;
            self.counter    =  0;
        end
        
        % resetting the model
        function reset(self)
            self.W_inh      = self.Cprob(meshgrid(self.Theta)...
                -meshgrid(self.Theta)',...
                self.a_i,self.c_i);
            self.Eta        = self.eta*ones(self.N);
            self.OD         = self.OD_0;
        end
        
        % setting orientation difference
        function set_OD(self,varargin)
            p               = inputParser;
            addOptional(p,'OD',[]);
            p.parse(varargin{:});
            od              = p.Results.OD;
            if ~isempty(od)
                self.OD         = od;
            else
                self.OD         = self.OD_0;
            end
        end
        
        % setting reference orientation (phi)
        function set_PHI(self,varargin)
            p               = inputParser;
            addOptional(p,'PHI',[]);
            p.parse(varargin{:});
            phi             = p.Results.PHI;
            if ~isempty(phi)
                self.Phi    = phi;
            else
                self.Phi    = self.Phi_0;
            end
        end
        
        % setting pp_correct
        function set_CORRECT(self,p)
           self.p_correct = p; 
        end
        
        % setting pp_correct
        function set_INCORRECT(self,p)
           self.p_incorrect = p; 
        end
        
        % getting just noticeable difference
        function JND = get_JND(self)
            JND             = self.mean_JND;
        end
        
        function [We,Wi] = get_weights(self)
            We              = self.W_exc;
            Wi              = self.W_inh;
        end
        
        function [v,r] = get_response(self)
            W               = self.W_exc-self.W_inh;
            V_ff            = self.J_ff*exp(...
                -((self.Adiff(self.Theta,self.Phi)).^2)...
                /(2*self.sigma_ff^2));
            [~,v]           = ode45(@(t,v)self.dV(v,V_ff,...
                W,self.alpha,self.tau),...
                [0 self.t_sim],self.V_0);
            v               = v(end,:)';
            r               = self.alpha*max(v,0);
        end
        
        % simulation of training session
        function session(self)
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
            W               = self.W_exc-self.W_inh;
            V_ff            = self.J_ff*exp(...
                -((self.Adiff(self.Theta,self.Phi)).^2)...
                /(2*self.sigma_ff^2));
            [~,v]           = ode45(@(t,v)self.dV(v,V_ff,...
                W,self.alpha,self.tau),...
                [0 self.t_sim],self.V_0);
            r               = self.alpha*max(v(end,:),0)';
            M_ref           = r*self.t_sim;
            Phi_probe       = self.Phi+((rand()>.5)*2-1)*self.OD;
            V_ff            = self.J_ff*exp(...
                -((self.Adiff(self.Theta,Phi_probe)).^2)...
                /(2*self.sigma_ff^2));
            [~,v]           = ode45(@(t,v)self.dV(v,V_ff,...
                W,self.alpha,self.tau),...
                [0 self.t_sim],self.V_0);
            r               = self.alpha*max(v(end,:),0)';
            M_probe         = r*self.t_sim;
            D               = abs(M_ref-M_probe)./...
                ((self.k*(M_ref+M_probe)).^.5);
            p               = .5*erfc(-D/(2^.5));
            p(isnan(p))     = .5;
            correct         = mean(p>rand(self.N,1))>=self.C;
            
            if ~correct
                if rand() < (1 - self.p_correct)
                    dW_inh          = self.eta*((1-self.W_inh/...
                        (self.J_rec*self.c_i*2^self.a_i)).^self.mu).*...
                        (r*r');
                    self.W_inh      = self.W_inh+dW_inh;
                end
                self.OD         = self.OD*1.2;
                self.counter    = 1;
            else
                if rand() < self.p_incorrect
                    dW_inh          = self.eta*((1-self.W_inh/...
                        (self.J_rec*self.c_i*2^self.a_i)).^self.mu).*...
                        (r*r');
                    self.W_inh      = self.W_inh+dW_inh;
                end
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