
function [outAll] = aplicaCOFAll(sinal,pt)

H=[ 1.0000    0.5633    0.1493    0.0424 0 0 0; ...
        0.4524    1.0000    0.5633    0.1493    0.0424 0 0; ...
        0.0172    0.4524    1.0000    0.5633    0.1493    0.0424 0; ...
        0    0.0172    0.4524    1.0000    0.5633    0.1493    0.0424; ...
        0 0    0.0172    0.4524    1.0000    0.5633    0.1493; ...
        0 0 0    0.0172    0.4524    1.0000    0.5633; ...
        0 0 0 0    0.0172    0.4524    1.0000];
n=7;
cent=4;

outDM=(sinal)*inv(H);
outCOF=zeros(size(outDM,1),1);
outAll=zeros(size(outDM,1),7);
for i=1:size(outDM,1)
    
   pu=zeros(1,n);
   pu(cent)=1;
   aux=[];
   det=0;
   for j=1:n
       if j==cent continue; end;
       if outDM(i,j)>pt
           pu(j)=1;
       end
   end
   pos = find(pu>0);

   D=[];
   count=0;
   for j=1:n    
      if pu(j)==1
         count=count+1;
         if (j==cent)
             ind=count;
         end
         D=[D;H(j,:)]; 
      end   
   end
   
   w=inv(D*D')*D;
   out=(sinal(i,:))*w';
   outAll(i,pos) = out; 
   outCOF(i)=out(ind);
   
end

end
