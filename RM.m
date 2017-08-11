classdef RM < handle % Mario Senden (10-08-2017)
    properties (Access = public)
        JND                     % just noticeable difference
        W                       % connection weight matrix
        Phi                     % stimulus orientation
    end
    properties (Access = private)
        % functions
        f_V                     % neuron dynamics       
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
        dt                      % time step in numerical integration
        tau                     % membrane time constant
        t_steps                 % number of simulation steps
        t_span                  % time span vector
        Theta                   % preferred orientation of each neuron
        V_0                     % initial values of membrane potential
        
    end
    methods
        % constructor
        function self = RM(varargin) 
            p = inputParser;
            addOptional(p,'N',512);
            addOptional(p,'sigma_ff',45);
            addOptional(p,'alpha',10);
            addOptional(p,'J_ff',.5);
            addOptional(p,'J_rec',1);
            addOptional(p,'a_e',2.2);
            addOptional(p,'a_i',1.4);
            addOptional(p,'c_e',.00118);
            addOptional(p,'c_i',.00165);
            addOptional(p,'k',1.47);
            addOptional(p,'C',.53);
            addOptional(p,'eta',15e-9);
            addOptional(p,'t_sim',.5);
            addOptional(p,'dt',2e-3);
            addOptional(p,'tau',15e-3);
            addOptional(p,'JND',7.5);
            
            p.parse(varargin{:});
            
            self.f_V = @(V_rec,V_ff,W,alpha,tau)...
            (-V_rec+V_ff+W*alpha*max(V_rec,0))/tau;
            self.Adiff       = @(A,B)...
            angle(exp((A-B)*pi/90*1i))*90/pi;
            self.Cprob       = @(In,a,c)...
            c*(cosd(2*In)+1).^a;
            
            self.N          = p.Results.N;
            self.sigma_ff   = p.Results.sigma_ff;
            self.alpha      = p.Results.alpha;
            self.J_ff       = p.Results.J_ff;
            self.J_rec      = p.Results.J_rec; 
            self.a_e        = p.Results.a_e;
            self.a_i        = p.Results.a_i;
            self.c_e        = p.Results.c_e;
            self.c_i        = p.Results.c_i; 
            self.k          = p.Results.k;
            self.C          = p.Results.C;
            self.eta        = p.Results.eta;
            self.dt         = p.Results.dt;
            self.tau        = p.Results.tau;
            self.t_sim      = p.Results.t_sim;
            
            
            self.t_steps    = floor(self.t_sim/self.dt)+1;  
            self.t_span     = linspace(0,self.t_sim,self.t_steps);
            self.V_0        = zeros(self.N,1);
            self.JND        = p.Results.JND;
            self.Theta      = linspace(-90,90,self.N)';
            self.Phi        = [];
            self.W          = (self.Cprob(meshgrid(self.Theta)...
                                -meshgrid(self.Theta)',...
                                self.a_e,self.c_e)...
                                -self.Cprob(meshgrid(self.Theta)...
                                -meshgrid(self.Theta)',...
                                self.a_i,self.c_i));
 
        end
        % network simulation
        function self = trial(self)
            V_ff            = self.J_ff*exp(...
                                -((self.Adiff(self.Theta,self.Phi)).^2)...
                                /(2*self.sigma_ff^2));
            [~,V_ref]       = ode45(@(t,V_ref)self.f_V(V_ref,V_ff,...
                                self.W,self.alpha,self.tau),...
                                self.t_span,self.V_0);
            M_ref           = self.alpha*max(V_ref(end,:),0)*self.t_sim;
                            
            Phi_probe       = self.Phi+((rand()>.5)*2-1)*self.JND;
            V_ff            = self.J_ff*exp(...
                                -((self.Adiff(self.Theta,Phi_probe)).^2)...
                                /(2*self.sigma_ff^2));
            [~,V_probe]     = ode45(@(t,V_probe)self.f_V(V_probe,V_ff,...
                                self.W,self.alpha,self.tau),...
                                self.t_span,self.V_0);           
            M_probe         = self.alpha*max(V_probe(end,:),0)*self.t_sim;
    
    
            self.JND=M_probe-M_ref; % the average JND per session needs to be tracked!
                                    % that differs of course from the exact
                                    % value to be achanged in the staircase
        end
        
    end
    
end