%This code plots and animates data from my N-body code.  I commented out
%the plots of Mercury, Venus, Uranus, and Neptune.
clear
clc

%This writes to an avi file that I converted to MP4 using EasyHtml5Video.
writerObj = VideoWriter('planet_alt.avi');
open(writerObj);

%Here I am loading the text file that I created in my N-body code and
%assigning each column of values in the file to a vector.
    load planet.txt;    %Sun
    x1=planet(:,1);
    y1=planet(:,2);
    z1=planet(:,3);
    
    x2=planet(:,4);     %Mercury
    y2=planet(:,5);
    z2=planet(:,6);
    
    x3=planet(:,7);     %Venus
    y3=planet(:,8);
    z3=planet(:,9);
    
    x4=planet(:,10);    %Earth
    y4=planet(:,11);
    z4=planet(:,12);
    
    x5=planet(:,13);    %Mars
    y5=planet(:,14);
    z5=planet(:,15);
    
    x6=planet(:,16);    %Jupiter
    y6=planet(:,17);
    z6=planet(:,18);
    
    x7=planet(:,19);    %Saturn
    y7=planet(:,20);
    z7=planet(:,21);
    
    x8=planet(:,22);    %Uranus
    y8=planet(:,23);
    z8=planet(:,24);
    
    x9=planet(:,25);    %Neptune
    y9=planet(:,26);
    z9=planet(:,27);
    
    
    %I am playing around with different values of step size dt.  Some
    %values of dt cause VideoWriter errors.
    dt = 5;
    total_time = size(x4,1); %This is the time it takes to form the total length of the vector x4.  All vectors are the same length.
    
    %This loop plots points along the position vectors x y z of the
    %planets.  It pauses along the way so that the entire vectors are not
    %instantly plotted.  This creates the animation.
    
    whitebg('k')
%     axis fill
    axis vis3d
%            axis([-10.0, 10.0, -10.0, 10.0, -10.0, 10.0])
    
    for i = 1:10000          
        

        %These are coordinates for each mass that is displayed in the
        %animation.  I converted floats to strings.
%         sun_x_coordinate = num2str(x1(i));
%         sun_y_coordinate = num2str(y1(i));
%         sun_z_coordinate = num2str(z1(i));
%         
%         earth_x_coordinate = num2str(x4(i));        
%         earth_y_coordinate = num2str(y4(i));
%         earth_z_coordinate = num2str(z4(i));
%         
%         mars_x_coordinate = num2str(x5(i));
%         mars_y_coordinate = num2str(y5(i));
%         mars_z_coordinate = num2str(z5(i));
%         
%         jupiter_x_coordinate = num2str(x6(i));
%         jupiter_y_coordinate = num2str(y6(i));
%         jupiter_z_coordinate = num2str(z6(i));
        
         
                
        %Here are the plots for one frame of the video.  For some reason,
        %putting all particles in the same call of plot3 allowed me to plot
        %just the current point of each particle with no trails.
         %plot3(x1(i),y1(i),z1(i),'.r',x2(i),y2(i),z2(i),'.b',x3(i),y3(i),z3(i),'.g',x4(i),y4(i),z4(i),'.y',x5(i),y5(i),z5(i),'.c',x6(i),y6(i),z6(i),'.m',x7(i),y7(i),z7(i),'.w',x8(i),y8(i),z8(i),'.b',x9(i),y9(i),z9(i),'.y','MarkerSize',5);        
        
    

         %I will not use these calls of plot3.
         plot3(x1(i),y1(i),z1(i),'.r','MarkerSize',20); %plots Sun
         hold on
         plot3(x2(i),y2(i),z2(i),'-b','MarkerSize',5);   %plots Mercury
         hold on
         plot3(x3(i),y3(i),z3(i),'-g','MarkerSize',5);  %plots Venus
         hold on
         plot3(x4(i),y4(i),z4(i),'.b','MarkerSize',5);  %plots Earth
         hold on
         plot3(x5(i),y5(i),z5(i),'.g','MarkerSize',5);  %plots Mars
         hold on
         plot3(x6(i),y6(i),z6(i),'.c','MarkerSize',5);   %plots Jupiter
         hold on
         plot3(x7(i),y7(i),z7(i),'-g','MarkerSize',5);  %plots Saturn
         hold on
         plot3(x8(i),y8(i),z8(i),'-c','MarkerSize',5);  %plots Uranus
         hold on
         plot3(x9(i),y9(i),z9(i),'-y','MarkerSize',5);   %plots Neptune 
        
         
         %I will use an "if" statement to change axis limits for some range
         %of i in the loop, to prevent the axis from randomly changing in the
         %video.
%           axis([-20.0, 20.0, -20.0, 20.0, -20.0, 20.0]) 
        
        %These are handles for the coordinates that are displayed.  I will probably only use 2 of them, for Sun and Earth, with all planets plotted.
%         h = text(1.0,1.2,4.5,['\color{red}Sun\color{white}  ', sun_x_coordinate, '   ', sun_y_coordinate, '   ', sun_z_coordinate]);
%         h2 = text(1.0,1.2,3.5,['\color{blue}Earth\color{white}  ', earth_x_coordinate, '   ', earth_y_coordinate, '   ', earth_z_coordinate]);
%         h3 = text(1.0,1.2,2.5,['\color{green}Mars\color{white}  ', mars_x_coordinate, '   ', mars_y_coordinate, '   ', mars_z_coordinate]);
%         h4 = text(1.0,1.2,1.5,['\color{lightBlue}Jupiter\color{white}  ', jupiter_x_coordinate, '   ', jupiter_y_coordinate, '   ', jupiter_z_coordinate]);

        movie_func(i) = getframe;  
        
        %camorbit needs it's own loop with getframe!
        %for j=1:0.1:30
        %if i==1
        %camorbit(0,1);
        %movie_func(i) = getframe;
        %end
        %end
        
        writeVideo(writerObj,movie_func(i));
        
        
        %This loop repeats each frame 10 times so the movie appears to play slowly.
         %for j = 1:10        
         %    writeVideo(writerObj,movie_func(i));
         %end
        
        %handles are deleted to prevent text being displayed over more
        %text.  I need to remember to delete 2 of these when I delete 2
        %handles above.
%         delete(h)       
%         delete(h2)
%         delete(h3)
%         delete(h4)
        
 
    end
     
   close(writerObj);     