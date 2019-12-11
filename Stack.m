classdef Stack < handle
    
    properties
        stuff
        top_index
        max_elems = 10000
    end
    
    methods
        function self = Stack()
            self.max_elems = 1000000;
            self.stuff = cell(1, self.max_elems);
            self.top_index = 1;
        end
        
        function push(self, obj)
            self.stuff(self.top_index) = {obj};
            self.top_index = self.top_index + 1;
            if self.top_index > self.max_elems
                fprintf('OVERFLOW %d\n', self.max_elems);
            end
        end
        
        function obj = pop(self)
            if self.size() > 0
                self.top_index = self.top_index - 1;
                obj = self.stuff{self.top_index};
            else
                error('Stack is empty');
            end
        end
        
        function clear(self)
            while self.size() > 0
                self.pop();
            end
        end
        
        function sz = size(self)
            sz = self.top_index - 1;
        end
    end
    
end

