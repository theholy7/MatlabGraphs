classdef GraphObj
    %GraphObj Creates a Graph Object
    %   A GraphObj recieves the number of Nodes, maximum number of Edges
    %   and Node distance, and creates a graph respecting those
    %   conditions. It outputs a Laplacian matrix of the graph - L - and a
    %   vector with the initial values of each node - initX.
    
    properties (SetAccess = private)
        %inputs
        numNodes; %number of nodes in Graph
        numMaxEdges; %maximum number of node connections
        nodeDistance; %maximum distance of connection (increases until all nodes are connected)
    end
    properties (Hidden = true)
        %outputs
        graphLaplacian;
        graphConnection;
        initX;
        nodePositions;
    end
    
    methods
        %Constructor
        function obj = GraphObj(n, maxE, nodeD)
            if nargin > 0
                tStart = tic;

                obj.numNodes = n;
                obj.numMaxEdges = maxE;
                obj.nodeDistance = nodeD;

                obj.initX = 10*rand(n,1);
                obj.nodePositions = 10*rand(2,n);

                obj = obj.createGraph();
                obj.graphConnection = obj.graphLaplacian ~=0;

                tElapsed = toc(tStart);
                sprintf('Graph in: %.1fs', tElapsed)
            end 
        end
    end
    
    methods (Access = private)
        %Create Graph Connections
        function obj = createGraph(obj)
            %init Laplacian Matrix
            obj.graphLaplacian = zeros(obj.numNodes, obj.numNodes);
            
       
            %init Identity Matrix
            I = eye(obj.numNodes);
            
            %flag to iterate while and iteration counter
            to_iterate = 1;
            iterationCounter = 0;
            
            while to_iterate == 1
                
                %iterate over all lines and collumns of the matrix L (all nodes)
                for k = 1:obj.numNodes
                    for l = 1:obj.numNodes
                        
                        %number of edges of each node
                        nodeKEdges = obj.graphLaplacian(k,k);
                        nodeLEdges = obj.graphLaplacian(l,l);
                        
                        %If node distance < r AND
                        %If node K edges less or equal to 5
                        %If node L edges less or equal to 5
                        %Ensure K-L == L-K node connection (upper part of L matrix)
                        %If node K and L are not connected
                        if (norm(obj.nodePositions(:,k)-obj.nodePositions(:,l)) < obj.nodeDistance) && ...
                                (nodeKEdges <= obj.numMaxEdges) && ...
                                (nodeLEdges <= obj.numMaxEdges) && ...
                                (abs(k-l) > 0) && ...
                                (obj.graphLaplacian(k,l) == 0)
                            
                            %create Edge K-L
                            e =  I(:,k)-I(:,l);
                            %create matrix with K-L and L-K and add to Edge Matrix L
                            obj.graphLaplacian = obj.graphLaplacian + e*e';
                            pause(0.01);
                        end;
                    end;
                end;
                
                s = svd(obj.graphLaplacian);
                
                if s(obj.numNodes-1) > 1e-3
                    to_iterate = 0;
                    break;
                end;
                
                %increase distance that nodes can connect
                obj.nodeDistance = obj.nodeDistance+0.1;
                
                iterationCounter = iterationCounter +1;
                if(iterationCounter > 1000)
                    break;
                end;
            end;
        end
        
    end
    
    methods (Access = public)
        %Get Node Positions
        function nodepos = getNodePositions(obj)
            nodepos = obj.nodePositions;
        end
        
        %Get Graph Laplacian Matrix
        function graphL = getGraphLaplacian(obj)
            graphL = obj.graphLaplacian;
        end
        
        %Get Graph Connections Matrix
        function graphCM = getGraphConnections(obj)
            graphCM = obj.graphConnection;
        end
        
        %Get initial node values
        function initx = getInitX(obj)
            initx = obj.initX;
        end
        
        %Plot Graph
        function plotGraph(obj)
            wb = waitbar(0, 'Please wait... Drawing');
            h = figure;
            set(h,'Visible','off')
            hold on;
            for l = 1:obj.numNodes
                for k = l:obj.numNodes
                    if obj.graphConnection(l, k) == 1
                        a = obj.nodePositions(:,k);
                        b = obj.nodePositions(:,l);
                        plot([a(1) b(1)],...
                            [a(2) b(2)],...
                            'LineStyle', '-', 'Color', 'blue',...
                            'LineWidth', 1,...
                            'Marker','o','MarkerFaceColor','red',...
                            'MarkerEdgeColor', 'red',...
                            'MarkerSize', 8 ...
                            );
                    end
                end
                waitbar(l/obj.numNodes)
            end
            hold off;
            delete(wb)
            set(h,'Visible','on')
            pause(0.01);
        end
        
        %Redefined plot Graph for new connect Matrix
        function plotCustomGraph(obj, connectMatrix)
            
            %assert that impossivle channels are kept at zero
            zeroentry = obj.graphConnection == 0;
            assert(isequal(obj.graphConnection(zeroentry), connectMatrix(zeroentry)),...
                'Connect Matrix uses impossible connections')
            
            %plot graph nodes and custom edges based on logical matrix
            wb = waitbar(0, 'Please wait... Drawing');
            h = figure;
            set(h,'Visible','off')
            hold on;
            for l = 1:obj.numNodes
                for k = l:obj.numNodes
                    if connectMatrix(l, k) == 1
                        a = obj.nodePositions(:,k);
                        b = obj.nodePositions(:,l);
                        plot([a(1) b(1)],...
                            [a(2) b(2)],...
                            'LineStyle', '-', 'Color', 'blue',...
                            'LineWidth', 1,...
                            'Marker','o','MarkerFaceColor','red',...
                            'MarkerEdgeColor', 'red',...
                            'MarkerSize', 8 ...
                            );
                    end
                end
                waitbar(l/obj.numNodes)
            end
            hold off;
            delete(wb)
            set(h,'Visible','on')
            pause(0.01);
        end
        
        %Plot convergence of Nodes
        function plotConvergence(obj, W)
            
            finalX = zeros(obj.numNodes, obj.numNodes);
            finalX(:,1) = obj.initX;
            
            for i = 2:3*obj.numNodes
                finalX(:,i)=W*finalX(:,i-1);
            end
            
            prompt = 'Plot how many nodes?\n';
            plots = input(prompt)
            if plots > 0 && plots < obj.numNodes
                
                figure;
                hold on;
                for i = 1:plots
                    plot(finalX(i,:))
                end
                hold off;
            end
            
        end
        
        %Mean value of Initial Nodes
        function avgNodes = initAverage(obj)
           avgNodes = mean(obj.initX); 
        end
        
    end
end

