function [] = isomatrix_quiver(A,varargin)
    
    p = inputParser;
    
    %% set up default values for optional parameters: ('Color' and 'Labels')
    color = [0,0,0];
    labels = {'','',''};
    
    % validation of user input color:
    vectorValidator = @(x) validateattributes(x,{'numeric'},{'size',[1,3]});
    addParameter(p,'Color',color,vectorValidator)
        
    % validation of user input labels:
    errorMsg1 = strcat('Labels error: please provide vector of size=',' ',num2str(3),').'); 
    errorMsg2 = 'Incorrect label formatting (must be cell-array).'; 
    labelLength = @(x) assert(length(x)==3,errorMsg1);
    labelType = @(x) assert(length(x)==3,errorMsg2);
    addParameter(p,'Labels',labels);
        
    % read in optional parameters    
    [nParams] = length(varargin);
    for param = 1:1:(nParams/2)
        ind = (param-1)*2 + 1;        
        if strcmp(varargin{ind}, 'Color')
            color=varargin{ind+1};
        elseif strcmp(varargin{ind}, 'Labels')
            labels=varargin{ind+1};
            labelLength(labels);
            labelType(labels);
        end
    end
    
    h = gcf;
    figure_number=h.Number;
    figure(figure_number); hold on;
    
    %% quiver mesh:
    gridlines = 20;
    step = 1/gridlines;
    step0 =step/2;
    x_grid = 0:step:(1+step); 
    y_grid = x_grid;
    [P,Q] = meshgrid(x_grid,y_grid);
    
    %% replicator equation dynamics:
    f1 = arrayfun(@(p, q) A.W1([p, q, 1 - p - q]), P, Q);
    f2 = arrayfun(@(p, q) A.W2([p, q, 1 - p - q]), P, Q);
    f3 = arrayfun(@(p, q) A.W3([p, q, 1 - p - q]), P, Q);
    phi = P.*f1 + Q.*f2 + (1 - P - Q).*f3;
    Z1 = P.*(f1-phi);
    Z2 = Q.*(f2-phi);
    Z3 = (1 - P - Q).*(f3-phi);

    %% coordinate transformation:
    X = (P./(tan(pi/3)) + (1-P-Q)./(sin(pi/3)))*sin(pi/3);
    Y = P * (1/2)*tan(60*pi/180);

    %% arrow length:
    len = 0.0375/gridlines*20;
    HL = 10;
    HW = 8;

    for i = 1:length(P)
        for j = 1:length(Q)                
            if ( P(i,j) < step0) || (Q(i,j) < step0) || ( (P(i,j)+Q(i,j)) > 1-step0)
                %% outside of simplex

            else

                %% calculate direction & magnitude
                V = Z1(i,j);
                U = (Z3(i,j)-Z2(i,j))*cos(pi/3); 
                mag = sqrt(V^2 + U^2);
                
                if mag > 0
                    %% set direction & magnitude
                    ah = annotation('arrow','headStyle','cback1','HeadLength',HL,'HeadWidth',HW);

                    set(ah,'parent',gca);
                    
                    set(ah,'position',[X(i,j), Y(i,j), len*U/mag, len*V/mag]);
                    set(ah,'Color',color);
                else
                   plot(X(i,j), Y(i,j),'.','MarkerSize',15,'Color',color);
                end
            end
        end
    end
    add_labels(labels);
end