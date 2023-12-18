#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <climits>
#include <queue>
#include <algorithm>
#include <map>

using namespace std;

class Graph {
private:
    int n; // Liczba wierzchołków
    vector<vector<int> > adj_matrix; // Macierz sąsiedztwa

public:
    Graph(int n) : n(n), adj_matrix(n, vector<int>(n, 0)) {}

    void add(int u, int v, int weight) 
    {
        adj_matrix[u][v] = weight;
        adj_matrix[v][u] = weight;
    }

    int get(int u, int v)
    {
        return adj_matrix[u][v];
    }

    void generate() 
    {
        srand(time(nullptr));

        for (int i = 1; i < n; ++i) 
        {
            int weight = rand() % 100 + 1;
            add(i - 1, i, weight);
        }

        int extra_edges = rand() % (n * (n - 1) / 2 - (n - 1)) + 1;

        for (int i = 0; i < extra_edges; ++i) {
            int u = rand() % n;
            int v = rand() % n;
            int weight = rand() % 100 + 1;
            if (u != v && adj_matrix[u][v] == 0) {
                add(u, v, weight);
            }
        }
    }

    void print() 
    {
        for (int i = 0; i < n; ++i) 
        {
            for (int j = 0; j < n; ++j) 
            {
                cout << adj_matrix[i][j] << " ";
                if(adj_matrix[i][j]<10)
                    cout<<" ";
            }
            cout << endl;
        }
    }

    int value()
    {
        int sum = 0;
        for(int i=0; i < n; i++)
            for(int j = 0; j < n; j++)
                sum += adj_matrix[i][j];
        return sum;
    }

    Graph copy()
    {
        Graph p(n);
        for(int i = 0; i < n; i++)
            for(int j = 0; j < n; j++)
                p.add(i,j,adj_matrix[i][j]);
        return p;
    }

    int shortestPath(int start, int end)
    {
        int sum = 0;
        vector<int> road(n, -1);
        int current = 0;
        road[current] = start;
        adj_matrix[start][start] = -1;
        while(road[current] != end)
        {
            int low = INT_MAX, next = -1;
            for(int i = 0; i < n; i++)
                if(adj_matrix[road[current]][i] < low && adj_matrix[road[current]][i] != 0 && adj_matrix[i][i] != -1)
                {
                    low = adj_matrix[road[current]][i];
                    next = i;
                }
            if(next == -1)
            {
                current--;
            }
            else
            {
                adj_matrix[next][next] = -1;
                road[++current] = next;
            }
        }
        for(int i = 0; i < n; i++)
            adj_matrix[i][i] = 0;
        for(int i = 0; i < n -1; i++)
            if(road[i+1] != -1)
                sum += adj_matrix[road[i]][road[i+1]];
            else
                break;
        return sum;
    }

    vector<int> primMST(Graph& z, int m, vector<bool> terminale, bool ulosowiony) 
    {
        vector<int> MST(m);
        int pom = 0;
        if(ulosowiony)
            while(!terminale[pom]) pom = rand() % n;
        else
            while(!terminale[pom]) pom++;
        MST[0] = pom;
        z.add(pom,pom,-1);
        int curent = 0;

        for (int i = 0; i < m; i++)
        {
            int naj = INT_MAX, y;
                for (int j = 0; j < n; j++)
                    if (z.get(j,j) != -1 && naj > z.get(MST[curent],j) && terminale[j])
                    {
                        naj = z.get(MST[curent],j);
                        y = j;
                    }
            MST[++curent] = y;
            z.add(y,y,-1);
        }

        return MST;
    }

    void convertMSTtoSTEINER(Graph& z, vector<int> MST, int m, bool ulosowiony)
    {
        vector<bool> used(m, false);
        for(int i = 0; i < m-1; i++)
        {
            int start, koniec;
            if(ulosowiony)
            {
                int pom = rand() % m;
                while(used[pom]) pom = rand() % m;
                used[pom] = true;
                start = MST[pom];
                if(pom == m-1)
                    koniec = MST[0];
                else
                    koniec = MST[pom+1];
            }
            else
            {
                start = MST[i];
                koniec = MST[i+1];
            }

            vector<int> road(n, -1);
            int current = 0;
            road[current] = start;
            adj_matrix[start][start] = -1;
            while(road[current] != koniec)
            {
                int low = INT_MAX, next = -1;
                for(int j = 0; j < n; j++)
                    if(adj_matrix[road[current]][j] < low && adj_matrix[road[current]][j] != 0 && adj_matrix[j][j] != -1)
                    {
                        low = adj_matrix[road[current]][j];
                        next = j;
                    }
                if(next == -1)
                {
                    current--;
                }
                else
                {
                    adj_matrix[next][next] = -1;
                    road[++current] = next;
                }
            }
            for(int j = 0; j < n; j++)
                adj_matrix[j][j] = 0;
            for(int j = 0; j < n-1; j++)
                if(road[j+1] != -1)
                    z.add(road[j],road[j+1],adj_matrix[road[j]][road[j+1]]);
                else
                    break;
        }
    }

    Graph zachlanny(vector<bool> terminale, int m, bool ulosowiony)
    {
        Graph z(n);
        for(int i = 0; i < n; i++)
            for(int j = i+1; j < n; j++)
                if(terminale[i] && terminale[j])
                    z.add(i,j,shortestPath(i,j));

        vector<int> MST = primMST(z, m, terminale, ulosowiony);
        
        for(int i = 0; i < n; i++)
            for(int j = 0; j < n; j++)
                z.add(i,j,0);
        
        convertMSTtoSTEINER(z, MST, m, ulosowiony);

        return z;
    }

    bool Pathexist(int start, int end)
    {
        bool exist = true;
        vector<int> road(n, -1);
        int current = 0;
        road[current] = start;
        adj_matrix[start][start] = -1;

        while(road[current] != end)
        {
            int low = INT_MAX, next = -1;
            for(int i = 0; i < n; i++)
                if(adj_matrix[road[current]][i] < low && adj_matrix[road[current]][i] != 0 && adj_matrix[i][i] != -1)
                {
                    low = adj_matrix[road[current]][i];
                    next = i;
                }
            if(next == -1)
            {
                current--;
            }
            else
            {
                adj_matrix[next][next] = -1;
                road[++current] = next;
            }
            if(current == -1)
                {exist = false; break;}
        }

        for(int i = 0; i < n; i++)
            adj_matrix[i][i] = 0;
        
        return exist;
    }

    bool steinerowskie(Graph& z, vector<bool> terminale)
    {
        for(int i = 0; i < n; i++)
            for(int j = i+1; j < n; j++)
                if (!z.Pathexist(i,j) and terminale[i] and terminale[j])
                    return false;
        return true;
    }

    Graph silowy(vector<bool> terminale)
    {
        Graph best(n), solution(n), rekurencja(n);
        int min = INT_MAX;

        for(int i = 0; i < n; i++)
            for(int j = 0; j < n; j++)
                best.add(i, j, adj_matrix[i][j]);
        solution = best.copy();
        
        for(int i = 0; i < n; i++)
            for(int j = 0; j < n; j++)
                if(best.get(i,j) != 0)
                {
                    best.add(i,j,0);
                    if(steinerowskie(best, terminale))
                        {
                            rekurencja = best.silowy(terminale);
                            int current = rekurencja.value();
                            if(min > current)
                            {
                                min = current;
                                solution = rekurencja.copy();
                            }
                        }
                    best.add(i,j,adj_matrix[i][j]);
                }

        return solution;
    }

    Graph mutacja_a(Graph& z, vector<bool> terminale)
    {
        Graph nowy = z.copy();
        int v, u;
        while(steinerowskie(nowy, terminale))
        {
            v = rand() % n, u = rand() % n;
            nowy.add(u,v,0);
        }
        nowy.add(u,v,z.get(u,v));
        return nowy;
    }

    Graph mutacja_b(Graph& z, vector<bool> terminale)
    {
        Graph nowy = z.copy();
        int v, u;
        while(steinerowskie(nowy, terminale))
        {
            v = rand() % n, u = rand() % n;
            nowy.add(u,v,z.get(u,v));
            v = rand() % n; u = rand() % n;
            nowy.add(u,v,0);
        }
        nowy.add(u,v,z.get(u,v));
        return nowy;
    }

    Graph mutacja_c(vector<bool> terminale, Graph& g)
    {
        Graph nowy = g.copy();
        
        int v, u;
        while(steinerowskie(nowy, terminale))
        {
            v = rand() % n, u = rand() % n;
            nowy.add(u,v,g.get(u,v));
            v = rand() % n; u = rand() % n;
            while(nowy.get(u,v) == 0) {v = rand() % n; u = rand() % n;}
            nowy.add(u,v,0);
        }
        nowy.add(u,v,g.get(u,v));
        return nowy;
    }

    bool rowne(Graph& a, Graph& b, int n)
    {
        for(int i = 0; i < n; i++)
            for(int j = 0; j < n; j++)
                if(a.get(i,j) != b.get(i,j))
                    return false;
        return true;
    }

    Graph Krzyzowanie(Graph& a, Graph& b, vector<bool> terminale)
    {
        Graph potomek(n);
        for(int i = 0; i < n; i++)
            for(int j = 0; j < n; j++)
                if(a.get(i,j) == b.get(i,j))
                    potomek.add(i,j,a.get(i,j));
        while(!steinerowskie(potomek, terminale))
            {
                int v = rand() % n, u = rand() % n;
                potomek.add(v,u,a.get(v,u));
            }
        return potomek;
    }

    Graph metaheurystyka(vector<bool> terminale, Graph& g, int powtorzenia)
    {
        Graph best(n);
        for(int i = 0; i < n; i++)
            for(int j = 0; j < n; j++)
                best.add(i, j, adj_matrix[i][j]);
        
        Graph a = mutacja_a(best, terminale), b = mutacja_a(best, terminale), c = mutacja_b(best, terminale), d = mutacja_b(best, terminale), e = mutacja_c(terminale, g);

        int value_best = best.value(), value_a = a.value(), value_b = b.value(), value_c = c.value(), value_d = d.value(), value_e = e.value();
        std::map<int, Graph> values = {
        {value_best, best},
        {value_a, a},
        {value_b, b},
        {value_c, c},
        {value_d, d},
        {value_e, e}};

        std::vector<std::pair<int, Graph>> sorted_values(values.begin(), values.end());
        std::sort(sorted_values.begin(), sorted_values.end(), 
            [](const std::pair<int, Graph>& a, const std::pair<int, Graph>& b) {
                return a.first > b.first;
            });

        Graph docelowy = Krzyzowanie(sorted_values[0].second, sorted_values[1].second, terminale);
        if(rowne(best,docelowy,n)) powtorzenia++;
        if(powtorzenia == 5) return best;
        return docelowy.metaheurystyka(terminale, g, powtorzenia);
    }
};

int main() {
    int n = 6, m = 4; // Liczba wierzchołków w grafie i liczba terminali
    Graph g(n);

    g.generate();

    cout << "Wygenerowny graf:" << endl;
    g.print();

    //losowanie terminali
    vector<bool> terminale(n, false);
    for(int i = 0; i < m; i++)
    {
        int v = rand() % n;
        if(!terminale[v])
            terminale[v] = true;
        else
            i--;
    }

    cout<<"lista terminali:"<<endl;
    for(int i = 0; i < n; i++)
        if(terminale[i])
            cout << i << " ";

    cout << endl << "Drzewo Steinera algorytmem zachlannym:" << endl;
    Graph z = g.zachlanny(terminale, m, false);
    z.print();
    
    cout << "Drzewo Steinera ulosowionym algorytmem zachlannym:" << endl;
    Graph u = g.zachlanny(terminale, m, true);
    u.print();

    cout << "Drzewo Steinera algorytmem pelnego przegladu:" << endl;
    Graph b = g.silowy(terminale);
    b.print();

    cout << "Drzewo Steinera algorytmem z zakresu mataheurystyki (algorytm genetyczny):" << endl;
    Graph a = g.metaheurystyka(terminale, g, 0);
    a.print();

    return 0;
}
