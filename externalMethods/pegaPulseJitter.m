% funcao para ler pulsehi_physics.dat e pegar o sinal deterministico do
% tile inlcuindo o jitter uniforme [-12,12]
% bernardo smp 04/10/10

function [pulsehi] = pegaPulseJitter()
    
    pulso = load('pulsehi_physics.txt');
    pulso = [zeros(150,2);pulso;zeros(150,2)];

    jitter = randi([-6,6],1,1);
%     jitter = randi([-25,25],1,1);
%     jitter = round(4.*randn(1,1)); %jitter gaussiano
%     jitter = 0;
    % deformacao no pulso
%     def=0.02.*randn(1,7);
    def=0.05.*randn(1,7);
        
    zero = find(pulso(:,1)<0);
    zero = zero(end)+1;
    
    pulsehi = [pulso(zero-150+jitter,2)+def(1) pulso(zero-100+jitter,2)+def(2) pulso(zero-50+jitter,2)+def(3) ...
                pulso(zero+jitter,2)+def(4) pulso(zero+50+jitter,2)+def(5) pulso(zero+100+jitter,2)+def(6) pulso(zero+150+jitter,2)+def(7)];

end