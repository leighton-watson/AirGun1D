classdef Schrodinger < scheme.Scheme
    properties
        m % Number of points in each direction, possibly a vector
        h % Grid spacing
        x % Grid
        order % Order accuracy for the approximation

        D % non-stabalized scheme operator
        H % Discrete norm
        M % Derivative norm
        alpha

        D2
        Hi
        e_l
        e_r
        d1_l
        d1_r
        gamm
    end

    methods
        % Solving SE in the form u_t = i*u_xx -i*V;
        function obj = Schrodinger(m,xlim,order,V)
            default_arg('V',0);

            [x, h] = util.get_grid(xlim{:},m);

            ops = sbp.Ordinary(m,h,order);

            obj.D2 = sparse(ops.derivatives.D2);
            obj.H =  sparse(ops.norms.H);
            obj.Hi = sparse(ops.norms.HI);
            obj.M =  sparse(ops.norms.M);
            obj.e_l = sparse(ops.boundary.e_1);
            obj.e_r = sparse(ops.boundary.e_m);
            obj.d1_l = sparse(ops.boundary.S_1);
            obj.d1_r = sparse(ops.boundary.S_m);


            if isa(V,'function_handle')
                V_vec = V(x);
            else
                V_vec = x*0 + V;
            end

            V_mat = spdiags(V_vec,0,m,m);

            obj.D = 1i * obj.D2 - 1i * V_mat;

            obj.m = m;
            obj.h = h;
            obj.order = order;

            obj.x = x;
        end


        % Closure functions return the opertors applied to the own doamin to close the boundary
        % Penalty functions return the opertors to force the solution. In the case of an interface it returns the operator applied to the other doamin.
        %       boundary            is a string specifying the boundary e.g. 'l','r' or 'e','w','n','s'.
        %       type                is a string specifying the type of boundary condition if there are several.
        %       data                is a function returning the data that should be applied at the boundary.
        %       neighbour_scheme    is an instance of Scheme that should be interfaced to.
        %       neighbour_boundary  is a string specifying which boundary to interface to.
        function [closure, penalty] = boundary_condition(obj,boundary,type,data)
            default_arg('type','dirichlet');
            default_arg('data',0);

            [e,d,s] = obj.get_boundary_ops(boundary);

            switch type
                % Dirichlet boundary condition
                case {'D','d','dirichlet'}
                    tau = s * 1i*d;
                    closure = obj.Hi*tau*e';

                    switch class(data)
                        case 'double'
                            penalty = -obj.Hi*tau*data;
                        case 'function_handle'
                            penalty = @(t)-obj.Hi*tau*data(t);
                        otherwise
                            error('Wierd data argument!')
                    end

                % Unknown, boundary condition
                otherwise
                    error('No such boundary condition: type = %s',type);
            end
        end

        function [closure, penalty] = interface(obj,boundary,neighbour_scheme,neighbour_boundary)
            % u denotes the solution in the own domain
            % v denotes the solution in the neighbour domain
            [e_u,d_u,s_u] = obj.get_boundary_ops(boundary);
            [e_v,d_v,s_v] = neighbour_scheme.get_boundary_ops(neighbour_boundary);

            a =  -s_u* 1/2 * 1i ;
            b =  a';

            tau = b*d_u;
            sig = -a*e_u;

            closure = obj.Hi * (tau*e_u' + sig*d_u');
            penalty = obj.Hi * (-tau*e_v' - sig*d_v');
        end

        % Ruturns the boundary ops and sign for the boundary specified by the string boundary.
        % The right boundary is considered the positive boundary
        function [e,d,s] = get_boundary_ops(obj,boundary)
            switch boundary
                case 'l'
                    e = obj.e_l;
                    d = obj.d1_l;
                    s = -1;
                case 'r'
                    e = obj.e_r;
                    d = obj.d1_r;
                    s = 1;
                otherwise
                    error('No such boundary: boundary = %s',boundary);
            end
        end

        function N = size(obj)
            N = obj.m;
        end

    end

    methods(Static)
        % Calculates the matrcis need for the inteface coupling between boundary bound_u of scheme schm_u
        % and bound_v of scheme schm_v.
        %   [uu, uv, vv, vu] = inteface_couplong(A,'r',B,'l')
        function [uu, uv, vv, vu] = interface_coupling(schm_u,bound_u,schm_v,bound_v)
            [uu,uv] = schm_u.interface(bound_u,schm_v,bound_v);
            [vv,vu] = schm_v.interface(bound_v,schm_u,bound_u);
        end
    end
end