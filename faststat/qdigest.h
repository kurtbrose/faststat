/*
Several percentile algorithms


*/

#ifdef _WIN32
#define inline __inline
#endif

typedef struct {
    short next;
    short padddd;
    int min;
    unsigned long long count;
} Qdigest_node;


typedef struct {
    unsigned long long n;
    int k;
    char log2_k;
    short free_head;
    short free_tail;
    short num_free;
    short new_head;
    short new_tail;
    short generations[32];
    Qdigest_node *nodes;
} Qdigest;


static inline short qdigest_alloc(Qdigest *q) {
    short next;
    next = q->free_head;
    q->free_head = q->nodes[q->free_head].next;
    q->num_free--;
    return next;
}

static inline void qdigest_free(Qdigest *q, short index) {
    q->nodes[q->free_tail].next = index;
    q->free_tail = index;
    q->num_free++;
}

//merge two qnode lists, maintaining ascending order
static inline short merge_qnode_lists(Qdigest *q, short a, short b) {
    Qdigest_node *nodes;
    short head, tail, freed;
    nodes = q->nodes;
    if(nodes[a].min < nodes[b].min) { // a is smaller, it should be head
        head = a;
        a = nodes[a].next;
    } else if(nodes[a].min > nodes[b].min) { // b is smaller, it should be head
        head = b;
        b = nodes[b].next;
    } else {
            head = a;
            nodes[a].count += nodes[b].count;
            a = nodes[a].next;
            freed = b;
            b = nodes[b].next;
            qdigest_free(q, freed);
        }
    tail = head;
    while(a && b) {
        if(nodes[a].min < nodes[b].min) {
            nodes[tail].next = a;
            a = nodes[a].next;
        } else if(nodes[a].min > nodes[b].min) {
            nodes[tail].next = b;
            b = nodes[b].next;
        } else {
            nodes[tail].next = a;
            nodes[a].count += nodes[b].count;
            a = nodes[a].next;
            freed = b;
            b = nodes[b].next;
            qdigest_free(q, freed);
        }
        tail = nodes[tail].next;
    }
    //append any leftover nodes to the end of the list
    //if there aren't any leftover nodes, set tail.next to -1
    nodes[tail].next = b == -1 ? a : b;
    return head;
}

#define SORTERS_LEN 16

//given an initially unsorted list of qnodes, sort them and merge duplicates
static inline short sort_and_compress(Qdigest *q, short head) {
    Qdigest_node *nodes;
    short unsorted_head, next;
    short sorters[SORTERS_LEN];  
    // sorters[n] is the pointer to a sorted list of length 2 ** n, or -1
    int cur_sorter;
    for(cur_sorter=0; cur_sorter < SORTERS_LEN; cur_sorter++) { sorters[cur_sorter] = -1; }
    nodes = q->nodes;
    unsorted_head = head;
    while(unsorted_head != -1) { // combine same length lists until unsorted input exhausted
        if(sorters[0] == -1) {
            sorters[0] = unsorted_head;
            unsorted_head = nodes[unsorted_head].next;
            nodes[sorters[0]].next = -1;
        }
        if(unsorted_head == -1) { break; }
        next = unsorted_head;
        unsorted_head = nodes[unsorted_head].next;
        nodes[next].next = -1;
        for(cur_sorter = 0; cur_sorter < SORTERS_LEN && sorters[cur_sorter] != -1; cur_sorter++) {
            next = merge_qnode_lists(q, next, sorters[cur_sorter]);
        }
        assert(cur_sorter < SORTERS_LEN);
        sorters[cur_sorter] = next;
    }
    // combine differently lengthed lists till there is only one
    // find the first valid linked list in sorters
    for(cur_sorter = 0; cur_sorter < SORTERS_LEN; cur_sorter++) {
        if(sorters[cur_sorter] != -1) {
            next = sorters[cur_sorter];
            break;
        }
    }
    // roll up all the elements in sorters
    for(; cur_sorter < SORTERS_LEN; cur_sorter++) {
        if(sorters[cur_sorter] == -1) { continue; }
        next = merge_qnode_lists(q, sorters[cur_sorter], next);
    }
    return next;
}


static inline void compress_generation(Qdigest *q, short cur_nodes, short parents, char generation) {
    unsigned int mask;
    short cur, parents_cur, sibling, parent, parents_prev, prev, next;
    unsigned long long count;
    unsigned long long min_count;
    Qdigest_node *nodes;

    nodes = q->nodes;
    min_count = q->n / q->k;
    cur = cur_nodes;
    parents_cur = parents;
    prev = 0;
    parents_prev = 0;
    mask = 0xFFFFFFFF - ((1 << (generation + 2)) - 1);

    while(cur) {
        // test if current node has a sibling
        if(nodes[cur].next && (nodes[cur].min & mask) == (nodes[nodes[cur].next].min & mask)) {
            sibling = nodes[cur].next;
        } else {
            sibling = 0;
        }
        // walk parents to the last one that could possibly be cur parent
        while(parents_cur && (int)(nodes[cur].min & mask) > nodes[parents_cur].min) {
            parents_prev = parents_cur;
            parents_cur = nodes[parents_cur].next;
        }
        // test if current node has a parent
        if(parents_cur && nodes[parents_cur].min == (nodes[cur].min & mask)) {
            parent = parents_cur;
        } else {
            parent = 0;
        }
        next = nodes[sibling ? sibling : cur].next;
        // QDigest count is count(cur) + count(sibling) + count(parent)
        count = nodes[cur].count;
        if(sibling) { count += nodes[sibling].count; }
        if(parent) { count += nodes[parent].count; }
        if(count < min_count) {
            // remove node from current generation list
            if(prev) {
                nodes[prev].next = next;
            } else {
                cur_nodes = next;
            }
            // throw away sibling if one exists
            if(sibling) { qdigest_free(q, sibling); }
            // if there is a parent, throw away cur as well
            if(parent) {
                nodes[parent].count = count;
                qdigest_free(q, cur);
            } else {
                nodes[cur].count = count;
                // reduce granularity to next generation up's
                nodes[cur].min &= mask;
                nodes[cur].next = parents_cur;
                // splice into parents list if possible
                if(parents_prev) {
                    nodes[parents_prev].next = cur;
                } else { // if comes before current parents head, replace
                    parents = cur;
                }
                parents_prev = cur;
            }
            // check if we still need to free more nodes
            if(q->k <= q->num_free) { break; }
        } else { // if nothing merged, advance the prev pointer
            prev = sibling ? sibling : cur;
        }
        cur = next;
    }
    q->generations[generation] = cur_nodes;
    q->generations[generation + 1] = parents;
}

static inline void compress(Qdigest *q) {
    char i;
    short new_nodes;
    new_nodes = sort_and_compress(q, q->new_head);
    q->new_head = q->new_tail = 0;
    q->generations[0] = merge_qnode_lists(
        q, q->generations[0], new_nodes);
    for(i=0; q->num_free < q->k && i < 31; i++) {
        if(q->generations[i] == 0) { continue; }
        compress_generation(q, q->generations[i], q->generations[i+1], i);
    }
}

// this is probably very bad
static union converter {
    int intval;
    float floatval;
    char bytesval[4];
};

void qdigest_update(Qdigest *q, float x) {
    short new_node;
    union converter convert;
    if(!q->num_free) {
        compress(q);
    }
    new_node = qdigest_alloc(q);
    convert.floatval = x;
    q->nodes[new_node].min = convert.intval;
    if(q->new_tail) {
        q->nodes[q->new_tail].next = new_node;
        q->new_tail = new_node;
    } else {
        q->new_head = q->new_tail = new_node;
    }
}
