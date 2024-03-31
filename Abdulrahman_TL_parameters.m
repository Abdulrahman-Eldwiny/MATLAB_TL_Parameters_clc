%%%% Name : Abdulrahman Mamduh Mostafa Eldwiny
%%%% Section: 4

clear all
clc
Mu0=4*pi*(10)^-7;

%choose dc line or ac line 
voltage_type= input('dc or ac  ');

if strcmp(voltage_type,'dc')
    %Calculation of DC_resistance
    p=input('please enter the specific resistance of the line ');
    L=input('please enter the length of the line ');
    A=input('please enter the area of the line ');
    R_dc=(p*L)/A;
    fprintf('the dc resistance of the line == %.10f \n',R_dc);
    fprintf('there is no capacitance or inductance in dc transmission line');
else
    
%%%calculation of AC_resistance
fprintf('-> to calculate the AC resistance please enter the following values \n');
power_loss=input('->please input the power loss value of the transmission line  '); 
current=input('->please enter the current value moving in the transmission line  ');
R_AC=power_loss/(current^2);
fprintf('the value of the resistance of the line = %f ohm \n',R_AC)


disp('please enter one of the types give below ');
disp('(make sure to write the full sentence of the type and between quotation marks)!!!!')
disp('1- solid one wire conductor');
disp('2- single phase two wire conductor');
disp('3- composite conductors');
disp('4- three phase symmetrical conductor');
disp('5- three phase unsymmetrical conductor');
disp('6- three phase double solid conductors');
disp('7- three phase double bundle conductors');
type= input('');





if strcmp(type,'solid one wire conductor')
    disp('   |            ');
    disp('   |            ');
    disp('   |            ');
    disp('   |    d     P ');
    disp(' a |< - - - > O ');
    disp('   |            ');
    disp('   |            ');
    disp('   |            ');
    r= input('please enter the radius of the one wire conductor  ');
    d= input('please enter the distance from the point P to the wire a  ');
    GMR=r*.7788;
    GMD=d;
    L_internal=Mu0/(8*pi);
    L_external=(Mu0/(8*pi)) *log(GMD/GMR);
    L_tot=L_internal+L_external;
    fprintf('the value of the internal inductance of the line = %.10f  H \n ',L_internal);
    fprintf('the value of the external inductance of the line = %.10f H \n',L_external);
    fprintf('the value of the total inductance of the line = %.10f H \n',L_tot);
    C=(2*pi*8.85*10^(-12))/log(GMD/r);
    disp('Capacitance between the Transmission lines: ');
    disp(C);
    
    
    
elseif (strcmp(type,'single phase two wire conductor'))
    
    disp(' |      | ');
    disp(' |      | ');
    disp(' |      | ');
    disp(' |   d  | ');
    disp('a|<- - >| b');
    disp(' |      | ');
    disp(' |      | ');
    disp(' |      | ');
    r= input('please enter the radius of the two similar wire conductor  ');
    GMD= input('please enter the distance between the two wire  d');
    GMR=r*.7788;
    
    L=(Mu0/pi)*log(GMD/GMR);
    fprintf('the value of the total inductance of the line = %.10f H \n',L);
    C=(2*pi*8.85*10^(-12))/log(GMD/r);
    disp('Capacitance between the Transmission lines:');
    disp(C);
%%%%%%%%%%%%%%%%%% COMPOSITE 



elseif strcmp(type,'composite conductors')
    number_of_groups=input('input 1 for 1 group and 2 for 2 group of conductors');
    %one group of conductors with a point (stranded)
        if (number_of_groups==1)
            
            disp(' 1 O             ')
            disp(' 2 O         p   ')
            disp(' 3 O        O    ')
            disp(' 4 O             ')
            disp(' n O             ')
            n=input('the number of strands n');
            m=1; %point
            r=input('please enter the radius of every strand');
            
            GMRL=(0.7788*r)^n;
            GMRC=(r)^n;
            GMD=1;
            
            for i=1:n
            
                for j=1:n-1
                    
                    fprintf('the distance between the %i strand and the next %i th in same group is  ',i,j)
                    distance_strands=input('');
                    GMRL=GMRL*distance_strands;
                    GMRC=GMRC*distance_strands;
                end
                    
                d=input('please enter the distance between the strand and the point');
                GMD=GMD*d;
            end
            
            GMDT=GMD^1/(n);
            GMRLT=(GMRL)^1/(n*n);
            GMRCT=(GMRC)^1/(n*n);
            L=2*10^(-7)*log(GMDT/GMRLT);
            fprintf('the value of GMR is %f',GMRT);
            fprintf('the value of the total inductance of the conductor = %f ',L);
            C=(2*pi*8.85*10^(-12))/log(GMDT/GMRCT);
            disp('Capacitance between the Transmission lines:');
            disp(C);
           
           %%%%%%%%%%%%%2 GROUPS
        else
                disp(' 1 O       1 O    ')
                disp(' 2 O       2 O    ')
                disp(' 3 O       3 O    ')
                disp(' 4 O       4 O    ')
                disp(' n O       m O    ')
               
                n=input('the number of strands of the first conductor');
                m=input('the number of strands of the second conductor');
                r1=input('please enter the radius of every strand in first conductor');
                r2=input('please enter the radius of every strand in second conductor');
 
                GMRL1=(r1*.7788)^n;
                GMRL2=(r2*.7788)^m;
                GMRC1=r1^n;
                GMRC2=r2^m;
                GMD1=1;
                GMD2=1;
                %%CALCULATION OF LA
                for i=1:n
                    
                    for j=1:n-1
                        fprintf('the distance between the %i strand and the next %i th in same group is  ',i,j)
                        distance_strands=input('');
                        GMRL1=GMRL1*distance_strands;
                        GMRC1=GMRC1*distance_strands;
                    end
                    
                    for k=1:m
                        fprintf('the distance between the %i strand in first group and the %i in second group  ',i,k)
                        distance_between_groups=input('');
                        GMD1=GMD1*distance_between_groups;
                    end
                    
                end
                GMRTL1=(GMRL1)^(1/(n*n));
                GMRTC1=(GMRC1)^(1/(n*n));
                GMDT1=GMD1^(1/(n*m));
                LA=2*10^(-7)*log(GMDT1/GMRTL1);
                
                fprintf('GMRTL1 %f \n',GMRTL1);
                fprintf('GMDT %f \n ',GMDT1);
                
                fprintf('the value of the inductance of the first group of conductors = %f \n ',LA);
                CA=(2*pi*8.85*10^(-12))/log(GMD1/GMRC1);
                disp('Capacitance between the Transmission lines: \n');
                disp(CA);
                
                %%% calculation of LB
                
                for i=1:m
                    for j=1:m-1
                        fprintf('the distance between the %i strand and the next %i th in same group  ',i,j);
                        distance_strands=input('');
                        GMRL2=GMR2*distance_strands;
                        GMRC2=GMR2*distance_strands;
                    end
                    for k=1:n
                        fprintf('the distance between the %i strand in first group and the %i in second group',i,k);
                        distance_between_groups=input('');
                        GMD2=GMD2*distance_between_groups;
                    end
                end
                GMRTL2=(GMRL2)^(1/(m*m));
                GMRTC2=(GMRC2)^(1/(m*m));

                GMDT2=GMD^(1/(n*m));
                LB=2*10^(-7)*log(GMDT2/GMRTL2);
                
                fprintf('the value of the inductance of the second group of conductors = %f ',LB);
                CB=(2*pi*8.85*10^(-12))/log(GMD2/GMRC2);
                disp('Capacitance between the Transmission lines:');
                disp(CB);
                
                Ltot=LA+LB;
                fprintf('the value of the total inductance between the two groups of conductors = %f ',Ltot);
                C=CA+CB;
                disp('Total Capacitance between the Transmission lines:');
                disp(C);
        end
         
        
        
%%%%%%%%%%%%%%%%%%%%% 3 phase unsymmetrical;

       elseif strcmp(type,'three phase unsymmetrical conductor')
           bundle_state=input('is every phase bundled or not (write yes or no) ');
           if strcmp(bundle_state,'no')
               
               %%unsymmetrical  not bundled
               
            disp('3.THREE PHASE SINGLE CIRCUIT (only for same diameters))');
            disp('      a 0 ')
            disp('      /   \')
            disp('     /     \')
            disp('    /       \')
            disp('b  0---------0 c')
            Dab=input('Enter Dab: ');
            Dbc=input('Enter Dbc: ');
            Dca=input('Enter Dca: ');
            d=input('Enter the diameter in meters: ');
            GMD=(Dab*Dbc*Dca)^(1/3);
            GMRL=(d/2)*exp(-1/4);
            GMRC=d/2;
            L=2*10^-7*log(GMD/GMRL);
            disp('Inductance in H/m: ');
            disp(L);
            C=(2*pi*8.85*10^(-12))/log(GMD/GMRC);
            disp('Capacitance between the Transimission lines in F/m: ');
            disp(C);
           else
               
               
               
               %% unsymmetrical bundled
            disp('3.THREE PHASE SINGLE CIRCUIT (only for same diameters))');
            disp('      a 0-0       ')
            disp('      /    \      ')
            disp('     /      \     ')
            disp('    /        \    ')
            disp(' b 0-0-------0-0 c')
            Dab=input('Enter Dab: ');
            Dbc=input('Enter Dbc: ');
            Dca=input('Enter Dca: ');
            d=input('Enter the diameter in meters: ');
            space=input('Enter the space between the bundles');
            GMD=(Dab*Dbc*Dca)^(1/3);
            GMRL=sqrt((d/2)*space*exp(-1/4));
            GMRC=sqrt((d/2)*space);
            L=2*10^-7*log(GMD/GMRL);
            disp('Inductance in H/m: ');
            disp(L);
            C=(2*pi*8.85*10^(-12))/log(GMD/GMRC);
            disp('Capacitance between the Transimission lines in F/m: ');
            disp(C);
           end
           
           
     %%%%%%%%%%%%%%%%%3phase symmetrical
     
        elseif strcmp(type,'three phase symmetrical conductor')
            bundle_state=input('is every phase bundled or not (write yes or no) ');
            if strcmp(bundle_state,'no')
                
                
                %%symmetrical not bundled
            disp('3.THREE PHASE SINGLE CIRCUIT (only for same diameters))');
            disp('     a 0')
            disp('      /  \')
            disp('     /    \')
            disp('    /      \')
            disp('   /        \')
            disp('b 0----------0 c')
            D=input('Enter one of the distances between two lines: ');
            d=input('Enter the diameter in meters: ');
            GMD=D;
            GMRL=(d/2)*exp(-1/4);
            GMRC=d/2;
            L=2*10^-7*log(GMD/GMRL);
            disp('Inductance in H/m: ');
            disp(L);
            C=(2*pi*8.85*10^(-12))/log(GMD/GMRC);
            disp('Capacitance between the Transmission lines:');
            disp(C);
            else
                
                
                %%%symmetrical bundled
            disp('3.THREE PHASE SINGLE CIRCUIT (only for same diameters))');
            disp('      a 0-0       ')
            disp('      /    \      ')
            disp('     /      \     ')
            disp('    /        \    ')
            disp(' b 0-0-------0-0 c')
            D=input('Enter one of the distances between two lines: ');
            d=input('Enter the diameter in meters: ');
            GMD=D;
            space=input('Enter the space between the bundles');
            GMRL=sqrt((d/2)*space*exp(-1/4));
            GMRC=sqrt((d/2)*space);
            L=2*10^-7*log(GMD/GMRL);
            disp('Inductance in H/m: ');
            disp(L);
            C=(2*pi*8.85*10^(-12))/log(GMD/GMRC);
            disp('Capacitance between the Transmission lines:');
            disp(C);
            end
            
            %%%%%%%%%%%%three phase double
            %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% three phase double --> not bundle conductors

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       elseif strcmp(type,'three phase double solid conductors')

                       fprintf('not necessary that y be larger than x\n');
        drawing = [
            
       '                           a1     x        c2       '
       ' Circuit Arrangement      .O-------------O          '
       ' -------------------- z1 . |             | .        '
       '                        .  |p            |  .       '
       '                      b1   |      y      |  b2      '
       '                        O-------------------O       '
       '                        .  |             | .        '
       '                         . |q            |.         '
       '                         c1|      z      | a2       '
       '                            O-----------O           '];
       disp(drawing)
           
        
        %Asking for the input parameters 
        r=input('enter the radius please');
        GMR=r*.7788;
        y1=input('\n\n ->Horizontal Distance of separation (x): ');
        y2=input('\n -> Horizontal Distance of separation (y): ');
        y3=input('\n -> Horizontal Distance of separation (z): ');
        y4=input('\n -> Vertical Distance of separation (p): ');
        y5=input('\n -> Vertical Distance of separation (q): ');

        if(y2>y1)
            y6=(y2-y1)/2;
        else
            y6=(y1-y2)/2;
        end  
            
        if(y3>y2)
            y7=(y3-y2)/2;
        else
            y7=(y2-y3)/2;
            
       
        end
        if(y1>y3)
            y8=(y1-y3)/2;
        else
            y8=(y3-y1)/2;
        end
        
        %Calculating individual distance between every conductor
        z1=(y4^2+(y6)^2)^(1/2); %a1b1
        z2=((y4+y5)^2+(y8)^2)^(1/2); % a1c1
        z3=((y4+y5)^2+(y3-y8)^2)^(1/2); %a1a2
        z4=((y2-y6)^2+(y4)^2)^(1/2); %a1b2
        z5=(y7^2+(y5)^2)^(1/2); % c1b1
        z6=((y2-y7)^2+(y4)^2)^(1/2); %a2b1

fprintf('\n\n * DISTANCE MEASUREMENTS :');
fprintf('\n\t|---------------------------|--------------|');
fprintf('\n\t| 1) D(a1b1,b1a1,b2c2,c2b2) |: %f m |',z1);
fprintf('\n\t| 2) D(a1c1,c1a1,a2c2,c2a2) |: %f m |',z2);
fprintf('\n\t| 3) D(a1a2,a2a1,c1c2,c2c1) |: %f m |',z3);
fprintf('\n\t| 4) D(a1b2,b2a1,b1c2,c2b1) |: %f m |',z4);
fprintf('\n\t| 5) D(c1b1,b1c1,b2a2,a2b2) |: %f  m  |',z5);
fprintf('\n\t| 6) D(a2b1,b1a2,b2c1,c1b2) |: %f m |',z6);
fprintf('\n\t| 7) D(a1c2,c2a1)           |: %f m |',y1);
fprintf('\n\t| 8) D(b1b2,b2b1)           |: %f m |',y2);
fprintf('\n\t| 1) D(c1a2,a2c1)           |: %f m |',y3);
fprintf('\n\t|___________________________|______________|'); 
n1=(z1*z4*z4*z5)^(1/4);
n2=(z1*z6*z6*z5)^(1/4);
n3=(z2*z2*y1*y3)^(1/4);
v=(n1*n2*n3)^(1/3);

fprintf('\n *GMD CALCULATION:');
fprintf('\n\t |-----------------------------------------------|');
fprintf('\n\t | 1)D(phase AB):                  %f m   |',n1);
fprintf('\n\t | 2)D(phase BC):                  %f m   |',n2);
fprintf('\n\t | 3)D(phase AC):                  %f m   |',n3);
fprintf('\n\t | Geometric Mean Distance (GMD) = <%f>m  |',v);
fprintf('\n\t |_______________________________________________|');
r1=GMR;
r2=r;
g1=(r1*z3)^(1/2);
g2=(r1*z3)^(1/2);
g3=(r1*y2)^(1/2);
h1=(r2*z3)^(1/2);
h2=(r2*z3)^(1/2);
h3=(r2*y2)^(1/2);
GMRL=(g1*g2*g3)^(1/3);       
GMRC=(h1*h2*h3)^(1/3);    
fprintf('\n *GMR CALCULATION:');
fprintf('\n\t|---------------------------------------------------------|');
fprintf('\n\t|Geometric Mean Radius for Inductor[GMRL]  =[%f] m  |',GMRL);              
fprintf('\n\t|Geometric Mean Radius for Capacitor[GMRC] =[%f] m  |',GMRC);       
fprintf('\n\t|_________________________________________________________|');

L=log(v/GMRL)*0.2;
fprintf('\n\n\t {INDUCTANCE (L)} =[%f] mH/km',L); 
C=0.0556/(log(v/GMRC));
fprintf('\n\t {CAPACITANCE (C)} =[%f] µF/km',C);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  three phase double --> bundle conductors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif strcmp(type,'three phase double bundle conductors ')
        frpintf('not necessary that y be larger than x');
        drawing = [
            
       '                           a1     x        c2       '
       ' Circuit Arrangement      .O-------------O          '
       ' -------------------- z1 . |             | .        '
       '                        .  |p            |  .       '
       '                      b1   |      y      |  b2      '
       '                        O-------------------O       '
       '                        .  |             | .        '
       '                         . |q            |.         '
       '                         c1|      z      | a2       '
       '                            O-----------O           '];
       disp(drawing)   
           
        
        %Asking for the input parameters 
        e=input('Its a 2 bundle conductor. Enter the Bundle spacing: ');
        r=input('enter the radius please');
        GMR=r*.7788;
        y1=input('\n\n ->Horizontal Distance of separation (x): ');
        y2=input('\n -> Horizontal Distance of separation (y): ');
        y3=input('\n -> Horizontal Distance of separation (z): ');
        y4=input('\n -> Vertical Distance of separation (p): ');
        y5=input('\n -> Vertical Distance of separation (q): ');

        if(y2>y1)
            y6=(y2-y1)/2;
        else
            y6=(y1-y2)/2;
        end  
            
        if(y3>y2)
            y7=(y3-y2)/2;
        else
            y7=(y2-y3)/2;
            
       
        end
        if(y1>y3)
            y8=(y1-y3)/2;
        else
            y8=(y3-y1)/2;
        end
        
        %Calculating individual distance between every conductor
        z1=(y4^2+(y6)^2)^(1/2); %a1b1
        z2=((y4+y5)^2+(y8)^2)^(1/2); % a1c1
        z3=((y4+y5)^2+(y3-y8)^2)^(1/2); %a1a2
        z4=((y2-y6)^2+(y4)^2)^(1/2); %a1b2
        z5=(y7^2+(y5)^2)^(1/2); % c1b1
        z6=((y2-y7)^2+(y4)^2)^(1/2); %a2b1

fprintf('\n\n * DISTANCE MEASUREMENTS :');
fprintf('\n\t|---------------------------|--------------|');
fprintf('\n\t| 1) D(a1b1,b1a1,b2c2,c2b2) |: %f m |',z1);
fprintf('\n\t| 2) D(a1c1,c1a1,a2c2,c2a2) |: %f m |',z2);
fprintf('\n\t| 3) D(a1a2,a2a1,c1c2,c2c1) |: %f m |',z3);
fprintf('\n\t| 4) D(a1b2,b2a1,b1c2,c2b1) |: %f m |',z4);
fprintf('\n\t| 5) D(c1b1,b1c1,b2a2,a2b2) |: %f  m  |',z5);
fprintf('\n\t| 6) D(a2b1,b1a2,b2c1,c1b2) |: %f m |',z6);
fprintf('\n\t| 7) D(a1c2,c2a1)           |: %f m |',y1);
fprintf('\n\t| 8) D(b1b2,b2b1)           |: %f m |',y2);
fprintf('\n\t| 1) D(c1a2,a2c1)           |: %f m |',y3);
fprintf('\n\t|___________________________|______________|'); 
n1=(z1*z4*z4*z5)^(1/4);
n2=(z1*z6*z6*z5)^(1/4);
n3=(z2*z2*y1*y3)^(1/4);
v=(n1*n2*n3)^(1/3);

fprintf('\n *GMD CALCULATION:');
fprintf('\n\t |-----------------------------------------------|');
fprintf('\n\t | 1)D(phase AB):                  %f m   |',n1);
fprintf('\n\t | 2)D(phase BC):                  %f m   |',n2);
fprintf('\n\t | 3)D(phase AC):                  %f m   |',n3);
fprintf('\n\t | Geometric Mean Distance (GMD) = <%f>m  |',v);
fprintf('\n\t |_______________________________________________|');
r1=(GMR*e)^(1/2);
r2=(r*e)^(1/2);
g1=(r1*z3)^(1/2);
g2=(r1*z3)^(1/2);
g3=(r1*y2)^(1/2);
h1=(r2*z3)^(1/2);
h2=(r2*z3)^(1/2);
h3=(r2*y2)^(1/2);
GMRL=(g1*g2*g3)^(1/3);       
GMRC=(h1*h2*h3)^(1/3);    
fprintf('\n *GMR CALCULATION:');
fprintf('\n\t|---------------------------------------------------------|');
fprintf('\n\t|Geometric Mean Radius for Inductor[GMRL]  =[%f] m  |',GMRL);              
fprintf('\n\t|Geometric Mean Radius for Capacitor[GMRC] =[%f] m  |',GMRC);       
fprintf('\n\t|_________________________________________________________|');

L=log(v/GMRL)*0.2;
fprintf('\n\n\t {INDUCTANCE (L)} =[%f] mH/km',L); 
C=0.0556/(log(v/GMRC));
fprintf('\n\t {CAPACITANCE (C)} =[%f] µF/km',C);

else
        fprintf('you add something else');
end
       
end
