close all

plotDir = 'D:\\WTF_EYE\\Plots';

locXcat1 = [locXYstruct(1).locX;locXYstruct(5).locX ];
locXcat2 = [locXYstruct(2).locX;locXYstruct(6).locX ];
locXcat3 = [locXYstruct(3).locX;locXYstruct(7).locX ];
locXcat4 = [locXYstruct(4).locX;locXYstruct(8).locX ];

locYcat1 = [locXYstruct(1).locY;locXYstruct(5).locY ];
locYcat2 = [locXYstruct(2).locY;locXYstruct(6).locY ];
locYcat3 = [locXYstruct(3).locY;locXYstruct(7).locY ];
locYcat4 = [locXYstruct(4).locY;locXYstruct(8).locY ];


h=figure('units','normalized','outerposition',[0 0 1 1]);

for p=1:4
    
    subplot(1,4,p);
    
    if      p==1; theseDataX = locXcat1; theseDataY = locYcat1; 
    elseif  p==2; theseDataX = locXcat2; theseDataY = locYcat2;
    elseif  p==3; theseDataX = locXcat3; theseDataY = locYcat3;
    elseif  p==4; theseDataX = locXcat4; theseDataY = locYcat4;
    end
    
    for j=1:size(theseDataX,1)
        for i=1:size(theseDataY,2)
            plot(theseDataX(j,i), theseDataY(j,i),'linestyle','none','marker','.','MarkerSize',10); hold on
            set(gca,'xlim',[0,1024],'ylim',[0,768])
        end
    end
    
    pbaspect([1,1,1])
    
     
end

saveas(h,[plotDir '/' 'Eye_Movements.jpeg'],'jpeg')
%saveas(h,[plotDir '/' 'Eye_Movements'],'epsc') not working
%saveas(h,[plotDir '/' 'Eye_Movements.pdf'],'pdf')





%     for j=1:size(locXYstruct(p).locX,1) %trials
%         for i=1:size(locXYstruct(p).locX,2) % 
%             plot(locXYstruct(p).locX(j,i),locXYstruct(p).locY(j,i),'linestyle','none','marker','.','MarkerSize',10); hold on
%             set(gca,'xlim',[0,1024],'ylim',[0,768])
%         end
%     end
%     