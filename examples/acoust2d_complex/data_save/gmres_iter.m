

% Data includes [ GMRES_count err(real)  err(imag)    N      E/proc    E_total,  total proc]


clear all    ;
format long e;

% flat exact: 
a= [
    67   0.103651203443800E-03     0.427159893200604E-03       3      16      16       1
   177   0.578346048829559E-07     0.237729332101810E-06       5      16      16       1
   283   0.218479678792960E-10     0.644367233182397E-10       7      16      16       1
   412   0.118840493001926E-10     0.645004188948306E-11       9      16      16       1
   546   0.144197986884365E-10     0.861367574489513E-11      11      16      16       1
   700   0.461870541812459E-10     0.172423568478575E-10      13      16      16       1
   851   0.338913341835223E-10     0.225021841437911E-10      15      16      16       1
];

% deform
b=[
   177  0.114784481371102E-02     0.777216345501652E-03       3      16      16       1
   331  0.540903995388531E-05     0.522454821016982E-05       5      16      16       1
   522  0.474404698858066E-07     0.395724022439481E-07       7      16      16       1
   733  0.495521068621940E-09     0.409915879018286E-09       9      16      16       1
   960  0.597085575426703E-10     0.392244015046117E-10      11      16      16       1
  1205  0.410685929708166E-10     0.687468970639316E-10      13      16      16       1
  1462  0.636367347706113E-10     0.724832416310051E-10      15      16      16       1
];

% double_exact
c=[
    77  0.221295016023415E-03     0.179347787844203E-03       3      16      16       1
   196  0.123478279462352E-06     0.107874318366896E-06       5      16      16       1
   322  0.358316709636597E-10     0.340371897333824E-10       7      16      16       1
   466  0.975203251485368E-11     0.975215654758221E-11       9      16      16       1
   618  0.119739773651872E-10     0.131514799051047E-10      11      16      16       1
   791  0.170197189675037E-10     0.188724019051234E-10      13      16      16       1
   959  0.300447999812548E-10     0.219487761299320E-10      15      16      16       1
];

% double_data  
d=[
   188  0.296414919516170E-02     0.230134002630314E-02       3      16      16       1
   362  0.123478279462352E-06     0.107874318366896E-06       5      16      16       1
   571  0.396148653880601E-07     0.319229496725804E-07       7      16      16       1
   801  0.662422339203772E-09     0.480078296560826E-09       9      16      16       1
  1044  0.119739773651872E-10     0.131514799051047E-10      11      16      16       1
  1304  0.103312497445884E-09     0.550658407760807E-10      13      16      16       1
  1577  0.114581399923708E-09     0.605563377220619E-10      15      16      16       1
];



for i=1:4;

   if (i==1) ; w=a; filename='gmres_flat_exact_alpha00'  ; end;
   if (i==2) ; w=b; filename='gmres_deform_data_alpha00' ; end;
   if (i==3) ; w=c; filename='gmres_double_exact_alpha00'; end;
   if (i==4) ; w=d; filename='gmres_double_data_alpha00' ; end;
   
   n = w(:,5).*(w(:,4)+1).^2;

   figure(i);set(gca,'fontsize',22);title(['GMRES: Iteration Count (dof= 2*n)']);
   figure(i);hold on;plot(n,w(:,1),'ks-','LineWidth',2,'MarkerSize',11,'MarkerEdgeColor','k','MarkerFaceColor','k')
   figure(i);xlabel('n=E(N+1)^2');ylabel('GMRES Iteration # ');axis([0 5000 0 2000]);
   figure(i);grid on;
   figure(i);print(filename,'-dpng' );
   figure(i);print(filename,'-depsc');

end


