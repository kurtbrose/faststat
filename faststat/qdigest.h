/*
Several percentile algorithms


*/

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
    q->free_head = q->nodes[free_head].next;
    q->num_free--;
    return next;
}

static inline void qdigest_free(Qdigest *q, short index) {
    q->nodes[free_tail].next = index;
    q->free_tail = index;
    q->num_free++;
}

static inline short nums2qnodes(Qdigest *q) {
    short ret;
    size_t i;
    Qdigest_node *nodes, head, cur;
    int cur_num;
    qsort();
    nodes = q->nodes;
    ret = _qdigest_alloc(q);
    head = nodes + ret;
    head->next = 0;
    head->min = q->input_buffer[0];
    head->count = 1;
    cur = head;
    for(i=1; i < q->k; i++) {
        cur_num = q->input_buffer[i];
        if(cur_num == cur->min) {
            cur->count++;
        } else {
            cur->next = _qdigest_alloc(q);
            cur = nodes + cur->next;
            cur->min = cur_num;
            cur->next = 0;
            cur->count = 1;
        }
    }
    return ret;
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
            nodes[a].count += nodes[b].count
            a = nodes[a].next;
            freed = b;
            b = nodes[b].next;
            _qdigest_free(q, freed);
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
            nodes[a].count += nodes[b].count
            a = nodes[a].next;
            freed = b;
            b = nodes[b].next;
            _qdigest_free(q, freed);
        }
        tail = nodes[tail].next;
    }
    //append any leftover nodes to the end of the list
    //if there aren't any leftover nodes, set tail.next to -1
    nodes[tail].next = b == -1 ? a : b;
    return head;
}

static inline short[2] compress_generation(Qdigest *q, short cur_nodes, short parents, char generation) {
    unsigned int mask;
    short cur, parents_cur, sibling, parent, parents_prev, prev;
    unsigned long long count;
    unsigned long long min_count;
    int target; // number to be freed
    Qdigest_node *nodes;

    target = q-> k - q->num_free;
    nodes = q->nodes;
    min_count = q->n / q->k;
    cur = cur_nodes;
    parents_cur = parents;
    prev = 0;
    parents_prev = 0;
    mask = 0xFFFFFFFF - ((1 << (generation + 2)) - 1);

    while(cur) {
        // test if current node has a sibling
        if(nodes[cur].next && nodes[cur].min & mask == nodes[nodes[cur].next] & mask) {
            sibling = nodes[cur].next;
        } else {
            sibling = 0;
        }
        // walk parents to the last one that could possibly be cur parent
        while(parents_cur && nodes[cur] & mask > nodes[parents_cur].min) {
            parents_prev = parents_cur;
            parents_cur = nodes[parents_cur].next;
        }
        // test if current node has a parent
        if(parents_cur && nodes[parents_cur].min == nodes[cur].min & mask) {
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
            if(sibling) { _qdigest_free(q, sibling); }
            // if there is a parent, throw away cur as well
            if(parent) {
                nodes[parent].count = count;
                _qdigest_free(q, cur);
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
            target--;
            if(!target) { break; }
        } else { // if nothing merged, advance the prev pointer
            prev = sibling ? sibling : cur;
        }
        cur = next;
    }

    return {cur_nodes, parents}
}

static inline void compress(QDigest *q) {
    char i;
    short new_heads[2];
    new_nodes = sort_and_compress(q, q->new_head);
    q->new_head = q->new_tail = 0;
    q->generations[0] = merge_qnode_lists(
        q, q->generations[0], new_nodes);
    for(i=0; q->num_free < q->k and i < 31; i++) {
        if(q->generations[gen_i] == 0) { continue; }
        new_heads = compress_generation(q, q->generations[i], q->generations[i+1], i);
        q->generations[i] = new_heads[0];
        q->generations[i + 1] = new_heads[1];
    }
}

// this is probably very bad
static union converter {
    int intval;
    float floatval;
    char bytesval[4];
}

void qdigest_update(Qdigest *q, float x) {
    short new_node;
    union converter convert;
    if(!q->num_free) {
        compress(q);
    }
    new_node = qdigest_alloc(q);
    convert.floatval = x;
    new_node.min = convert.intval;
    if(q->new_tail) {
        q->nodes[q->new_tail].next = new_node;
        q->new_tail = new_node;
    } else {
        q->new_head = q->new_tail = new_node;
    }
}
