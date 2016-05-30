addpath('/Users/onat/Documents/Code/Matlab/palamedes1_8_0/Palamedes/')
PF                  = @PAL_Weibull;
options             = PAL_minimize('options');
options.MaxIter     = 10.^6;
options.MaxFunEvals = 10.^6;
options.Display     = 'On';
options.TolX        = 10.^-4;
options.TolFun      = 10.^-4;
%%
sigmas   = linspace(10,90,10);%subject's precision of face representation in degrees
lambdas  = linspace(45,135,10);%subject's threshold
stims    = linspace(0,250,200);%location of stimuli
data     = NaN(length(stims),length(sigmas),length(lambdas));
params   = NaN(4,length(sigmas),length(lambdas));
lc       = 0;
%%
for lambda = lambdas(:)'
    lambda
    lc = lc + 1;sc = 0;
    for sig = sigmas(:)'
        sc = sc+1;cc = 0;        
        data(:,sc,lc) = 1-normcdf(lambda,stims,sig);
        %
        searchGrid.alpha  = 15;
        searchGrid.beta   = 2;%10.^pmf.beta(chain);
        searchGrid.gamma  = data(1,sc,lc);%pmf.gamma(chain);
        searchGrid.lambda = 1-data(end,sc,lc);%pmf.lambda(chain);
        params0 = [ searchGrid.alpha  searchGrid.beta searchGrid.gamma searchGrid.lambda];
        %        
        funny                                  = @(params) sum( (data(:,sc,lc)' - PF(params,stims)).^2);
        options                                = optimset('Display','off','maxfunevals',10000,'tolX',10^-12,'tolfun',10^-12,'MaxIter',10000,'Algorithm','interior-point');
        [o.params1, o.Likelihood, o.ExitFlag]  = fmincon(funny, params0, [],[],[],[],[-Inf -Inf 0 0],[Inf Inf 1 1],[],options);
        c = [sc./length(sigmas) 0 1-sc./length(sigmas)]
        figure(10);
        plot(data(:,sc,lc),'color',c);hold on;
        plot(PF(o.params1,stims),'o','color',c);drawnow;
        params(:,sc,lc) = o.params1;                
    end
    hold off;
    figure(20);
    for n = 1:3
        subplot(1,3,n);hold on;plot(params(n,:,lc),'o-','color',rand(1,3)); 
    end
    drawnow    
%     pause;
end
