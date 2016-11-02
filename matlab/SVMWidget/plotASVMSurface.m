function [tx,ty,tz, hx]= plotASVMSurface(learned_svm, input_data, face_alpha, plot_sv, plot_normals)

num_contours = 10;
N = size(learned_svm.Sva,1);

if(isempty(input_data) && plot_normals)
    warning('Cannot plot normals without a data set');
end

cla; hold on;
tz=[];
if(N==2)
    if(~isempty(input_data))
        hh=plot(input_data(1,:)',input_data(2,:)','k.');
        axis equal
        axis tight
        lims=[get(gca,'XLim'), get(gca,'YLim')];
        lims = lims - [(lims(2)-lims(1))*0.2 -(lims(2)-lims(1))*0.2 (lims(4)-lims(3))*0.2 -(lims(4)-lims(3))*0.2];
        delete(hh);
    else
        lims=[-1 1 -1 1];
    end
    [tx,ty]=meshgrid(linspace(lims(1), lims(2), 40),linspace(lims(3), lims(4), 40));
%     hx=[];
%     maxh=-inf;
%     minh=inf;
%     for i=1:size(tx,1)
%         for j=1:size(tx,2)
%             pt=[tx(i,j);ty(i,j)];
%             hx(i,j) = calculateClassifier(learned_svm, pt);
%             
%             if(hx(i,j) > maxh)
%                 maxh=hx(i,j);
%             end
%             if(hx(i,j) < minh)
%                 minh=hx(i,j);
%             end
%         end
%     end

        hx=[];
        pts=[];
        for i=1:size(tx,1)
            for j=1:size(tx,2)
                pts=[pts,[tx(i,j);ty(i,j)]];
            end
        end
        hxfull = mx_calculate_classifier(learned_svm, pts);
        maxh = max(hxfull);
        minh=min(hxfull);
        count=1;
        for i=1:size(tx,1)
            for j=1:size(tx,2)
                hx(i,j) = hxfull(count);
                count=count+1;
            end
        end
        
    
    pcolor(tx,ty,hx);
    contour(tx, ty, hx, linspace(minh, maxh, num_contours),'k-');
    contour(tx, ty, hx, 1e-10,'k--', 'Linewidth',2);
    if(~isempty(input_data))
        plot(input_data(1,:)',input_data(2,:)','k.');
    end
    
    if(plot_sv)
        for i=1:size(learned_svm.Svb,2)
            if(learned_svm.Svb(3,i) == 0)
                clr = 'g';
            else if (learned_svm.Svb(4,i) == 0)
                    clr = 'r';
                else
                    clr = 'k';
                end
            end
            line([learned_svm.Svb(1,i);learned_svm.Svb(1,i) + 0.1*learned_svm.Svb(3,i)] ,[learned_svm.Svb(2,i);learned_svm.Svb(2,i) + 0.1*learned_svm.Svb(4,i)],'Color',clr, 'Linewidth',2);
        end
        plot(learned_svm.Sva(1,:),learned_svm.Sva(2,:)','ko','Linewidth',2,'Markersize',15);
    end
    axis tight
    colormap autumn
    shading interp
else if(N==3)
        meshsize = 20;
        if(~isempty(input_data))
            hh=plot3(input_data(1,:)',input_data(2,:)',input_data(3,:)','k.');
            axis equal
            axis tight
            lims=[get(gca,'XLim'), get(gca,'YLim'), get(gca,'ZLim')];
            lims = lims - [(lims(2)-lims(1))*0.4 -(lims(2)-lims(1))*0.4 (lims(4)-lims(3))*0.4 -(lims(4)-lims(3))*0.4 (lims(6)-lims(5))*0.4 -(lims(6)-lims(5))*0.4];
            delete(hh);
        else
            lims = [-1 1 -1 1 -1 1];
        end
        [tx,ty,tz]=meshgrid(linspace(lims(1), lims(2), meshsize),linspace(lims(3), lims(4), meshsize),linspace(lims(5), lims(6), meshsize));
        
        hx=[];
        pts=[];
        for i=1:size(tx,1)
            for j=1:size(tx,2)
                for k=1:size(tx,3)
                    pts=[pts,[tx(i,j,k);ty(i,j,k);tz(i,j,k)]];
                end
            end
        end
        hxfull = mx_calculate_classifier(learned_svm, pts);
        maxh = max(hxfull);
        count=1;
        for i=1:size(tx,1)
            for j=1:size(tx,2)
                for k=1:size(tx,3)
                    hx(i,j,k) = hxfull(count);
                    count=count+1;
                end
            end
        end
        
        
        p = patch(isosurface(tx, ty, tz, hx,0));
        set(p,'FaceColor',[ 0.49; .49; .49],'EdgeColor','none','FaceAlpha',face_alpha);
        
%                 p = patch(isosurface(tx, ty, tz, hx,0.2));
%         set(p,'FaceColor',[ 0.49; .49; .49],'EdgeColor','none','FaceAlpha',face_alpha);
        
        if(~isempty(input_data))
            plot3(input_data(1,:)',input_data(2,:)',input_data(3,:)','k.', 'Markersize',12);
        end
        
        if(plot_sv)
            for i=1:size(learned_svm.Svb,2)
                if(learned_svm.Svb(4,i) > 0 && learned_svm.Svb(5,i) == 0 && learned_svm.Svb(6,i) == 0)
                    clr = 'r';
                else if (learned_svm.Svb(5,i) > 0  && learned_svm.Svb(4,i) == 0 && learned_svm.Svb(6,i) == 0)
                        clr = 'g';
                    else if (learned_svm.Svb(6,i) > 0  && learned_svm.Svb(5,i) == 0 && learned_svm.Svb(4,i) == 0)
                            clr='b';
                        else 
                            clr = 'k';
                        end
                    end
                end
                line([learned_svm.Svb(1,i);learned_svm.Svb(1,i) + 0.1*learned_svm.Svb(4,i)] ,...
                    [learned_svm.Svb(2,i);learned_svm.Svb(2,i) + 0.1*learned_svm.Svb(5,i)],...
                    [learned_svm.Svb(3,i);learned_svm.Svb(3,i) + 0.1*learned_svm.Svb(6,i)],'Linewidth',2,'Color',clr );
            end
            plot3(learned_svm.Sva(1,:),learned_svm.Sva(2,:)',learned_svm.Sva(3,:)','ko','Linewidth',2,'Markersize',15);
        end
        if(~isempty(input_data) && plot_normals)
            dhxfull = mx_calculate_classifier_derivative(learned_svm, input_data(1:N,:));
            for i=1:size(dhxfull,2)
                pt=input_data(1:N,i);
                tmpnrm = dhxfull(:,i)/norm(dhxfull(:,i));
                pt2 = pt + tmpnrm*0.2;
                plot3([pt(1); pt2(1)],[pt(2); pt2(2)],[pt(3); pt2(3)],'k-','Linewidth',2);
            end
        end
        
        axis tight
        grid on
        camlight; lighting phong;
    end
end